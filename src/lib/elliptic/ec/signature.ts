import { BN } from 'bn.js'
import { utils } from '../utils'

const assert = utils.assert

export function Signature(options, enc) {
	if (options instanceof Signature) return options

	if (this._importDER(options, enc)) return

	assert(options.r && options.s, 'Signature without r or s')
	this.r = new BN(options.r, 16)
	this.s = new BN(options.s, 16)
	if (options.recoveryParam === undefined) this.recoveryParam = null
	else this.recoveryParam = options.recoveryParam
}

function Position() {
	this.place = 0
}

function getLength(buf, p) {
	const initial = buf[p.place++]
	if (!(initial & 0x80)) {
		return initial
	}
	const octetLen = initial & 0xf

	// Indefinite length or overflow
	if (octetLen === 0 || octetLen > 4) {
		return false
	}

	let val = 0
	let i
	let off
	for (i = 0, off = p.place; i < octetLen; i++, off++) {
		val <<= 8
		val |= buf[off]
		val >>>= 0
	}

	// Leading zeroes
	if (val <= 0x7f) {
		return false
	}

	p.place = off
	return val
}

function rmPadding(buf) {
	let i = 0
	const len = buf.length - 1
	while (!buf[i] && !(buf[i + 1] & 0x80) && i < len) {
		i++
	}
	if (i === 0) {
		return buf
	}
	return buf.slice(i)
}

Signature.prototype._importDER = function _importDER(data, enc) {
	data = utils.toArray(data, enc)
	const p = new Position()
	if (data[p.place++] !== 0x30) {
		return false
	}
	const len = getLength(data, p)
	if (len === false) {
		return false
	}
	if (len + p.place !== data.length) {
		return false
	}
	if (data[p.place++] !== 0x02) {
		return false
	}
	const rlen = getLength(data, p)
	if (rlen === false) {
		return false
	}
	let r = data.slice(p.place, rlen + p.place)
	p.place += rlen
	if (data[p.place++] !== 0x02) {
		return false
	}
	const slen = getLength(data, p)
	if (slen === false) {
		return false
	}
	if (data.length !== slen + p.place) {
		return false
	}
	let s = data.slice(p.place, slen + p.place)
	if (r[0] === 0) {
		if (r[1] & 0x80) {
			r = r.slice(1)
		} else {
			// Leading zeroes
			return false
		}
	}
	if (s[0] === 0) {
		if (s[1] & 0x80) {
			s = s.slice(1)
		} else {
			// Leading zeroes
			return false
		}
	}

	this.r = new BN(r)
	this.s = new BN(s)
	this.recoveryParam = null

	return true
}

function constructLength(arr, len) {
	if (len < 0x80) {
		arr.push(len)
		return
	}
	let octets = 1 + ((Math.log(len) / Math.LN2) >>> 3)
	arr.push(octets | 0x80)
	while (--octets) {
		arr.push((len >>> (octets << 3)) & 0xff)
	}
	arr.push(len)
}

Signature.prototype.toDER = function toDER(enc) {
	let r = this.r.toArray()
	let s = this.s.toArray()

	// Pad values
	if (r[0] & 0x80) r = [0].concat(r)
	// Pad values
	if (s[0] & 0x80) s = [0].concat(s)

	r = rmPadding(r)
	s = rmPadding(s)

	while (!s[0] && !(s[1] & 0x80)) {
		s = s.slice(1)
	}
	let arr = [0x02]
	constructLength(arr, r.length)
	arr = arr.concat(r)
	arr.push(0x02)
	constructLength(arr, s.length)
	const backHalf = arr.concat(s)
	let res = [0x30]
	constructLength(res, backHalf.length)
	res = res.concat(backHalf)
	return utils.encode(res, enc)
}
