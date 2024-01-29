import { BN } from 'bn.js'
import { utils } from '../utils'
import { rand } from 'brorand'
import { KeyPair } from './key'
import { Signature } from './signature'
import { curves } from '../curves'
import { HmacDRBG } from 'hmac-drbg'

const assert = utils.assert

export class EC {
	curve: any
	n: any
	nh: any
	g: any
	hash: any

	constructor(options) {
		// Shortcut `elliptic.ec(curve-name)`
		if (typeof options === 'string') {
			assert(Object.prototype.hasOwnProperty.call(curves, options), 'Unknown curve ' + options)

			options = curves[options]
		}

		// Shortcut for `elliptic.ec(elliptic.curves.curveName)`
		if (options instanceof curves.PresetCurve) options = { curve: options }

		this.curve = options.curve.curve
		this.n = this.curve.n
		this.nh = this.n.ushrn(1)
		this.g = this.curve.g

		// Point on curve
		this.g = options.curve.g
		this.g.precompute(options.curve.n.bitLength() + 1)

		// Hash for function for DRBG
		this.hash = options.hash || options.curve.hash
	}

	keyPair(options) {
		return new KeyPair(this, options)
	}

	keyFromPrivate(priv, enc) {
		return KeyPair.fromPrivate(this, priv, enc)
	}

	keyFromPublic(pub, enc) {
		return KeyPair.fromPublic(this, pub, enc)
	}

	genKeyPair(options) {
		if (!options) options = {}

		// Instantiate Hmac_DRBG
		const drbg = new HmacDRBG({
			hash: this.hash,
			pers: options.pers,
			persEnc: options.persEnc || 'utf8',
			entropy: options.entropy || rand(this.hash.hmacStrength),
			entropyEnc: (options.entropy && options.entropyEnc) || 'utf8',
			nonce: this.n.toArray(),
		})

		const bytes = this.n.byteLength()
		const ns2 = this.n.sub(new BN(2))
		for (;;) {
			const priv = new BN(drbg.generate(bytes))
			if (priv.cmp(ns2) > 0) continue

			priv.iaddn(1)
			return this.keyFromPrivate(priv, undefined)
		}
	}

	_truncateToN(msg, truncOnly) {
		const delta = msg.byteLength() * 8 - this.n.bitLength()
		if (delta > 0) msg = msg.ushrn(delta)
		if (!truncOnly && msg.cmp(this.n) >= 0) return msg.sub(this.n)
		else return msg
	}

	sign(msg, key, enc, options) {
		if (typeof enc === 'object') {
			options = enc
			enc = null
		}
		if (!options) options = {}

		key = this.keyFromPrivate(key, enc)
		msg = this._truncateToN(new BN(msg, 16), undefined)

		// Zero-extend key to provide enough entropy
		const bytes = this.n.byteLength()
		const bkey = key.getPrivate().toArray('be', bytes)

		// Zero-extend nonce to have the same byte size as N
		const nonce = msg.toArray('be', bytes)

		// Instantiate Hmac_DRBG
		const drbg = new HmacDRBG({
			hash: this.hash,
			entropy: bkey,
			nonce: nonce,
			pers: options.pers,
			persEnc: options.persEnc || 'utf8',
		})

		// Number of bytes to generate
		const ns1 = this.n.sub(new BN(1))

		for (let iter = 0; ; iter++) {
			let k = options.k ? options.k(iter) : new BN(drbg.generate(this.n.byteLength()))
			k = this._truncateToN(k, true)
			if (k.cmpn(1) <= 0 || k.cmp(ns1) >= 0) continue

			const kp = this.g.mul(k)
			if (kp.isInfinity()) continue

			const kpX = kp.getX()
			const r = kpX.umod(this.n)
			if (r.cmpn(0) === 0) continue

			let s = k.invm(this.n).mul(r.mul(key.getPrivate()).iadd(msg))
			s = s.umod(this.n)
			if (s.cmpn(0) === 0) continue

			let recoveryParam = (kp.getY().isOdd() ? 1 : 0) | (kpX.cmp(r) !== 0 ? 2 : 0)

			// Use complement of `s`, if it is > `n / 2`
			if (options.canonical && s.cmp(this.nh) > 0) {
				s = this.n.sub(s)
				recoveryParam ^= 1
			}

			return Signature({ r: r, s: s, recoveryParam: recoveryParam }, undefined)
		}
	}

	verify = function verify(msg, signature, key, enc) {
		msg = this._truncateToN(new BN(msg, 16))
		key = this.keyFromPublic(key, enc)
		signature = Signature(signature, 'hex')

		// Perform primitive values validation
		const r = signature.r
		const s = signature.s
		if (r.cmpn(1) < 0 || r.cmp(this.n) >= 0) return false
		if (s.cmpn(1) < 0 || s.cmp(this.n) >= 0) return false

		// Validate signature
		const sinv = s.invm(this.n)
		const u1 = sinv.mul(msg).umod(this.n)
		const u2 = sinv.mul(r).umod(this.n)
		let p

		if (!this.curve._maxwellTrick) {
			p = this.g.mulAdd(u1, key.getPublic(), u2)
			if (p.isInfinity()) return false

			return p.getX().umod(this.n).cmp(r) === 0
		}

		// NOTE: Greg Maxwell's trick, inspired by:
		// https://git.io/vad3K

		p = this.g.jmulAdd(u1, key.getPublic(), u2)
		if (p.isInfinity()) return false

		// Compare `p.x` of Jacobian point with `r`,
		// this will do `p.x == r * p.z^2` instead of multiplying `p.x` by the
		// inverse of `p.z^2`
		return p.eqXToP(r)
	}

	recoverPubKey(msg, signature, j, enc) {
		assert((3 & j) === j, 'The recovery param is more than two bits')
		signature = Signature(signature, enc)

		const n = this.n
		const e = new BN(msg)
		let r = signature.r
		const s = signature.s

		// A set LSB signifies that the y-coordinate is odd
		const isYOdd = j & 1
		const isSecondKey = j >> 1
		if (r.cmp(this.curve.p.umod(this.curve.n)) >= 0 && isSecondKey)
			throw new Error('Unable to find sencond key candinate')

		// 1.1. Let x = r + jn.
		if (isSecondKey) r = this.curve.pointFromX(r.add(this.curve.n), isYOdd)
		else r = this.curve.pointFromX(r, isYOdd)

		const rInv = signature.r.invm(n)
		const s1 = n.sub(e).mul(rInv).umod(n)
		const s2 = s.mul(rInv).umod(n)

		// 1.6.1 Compute Q = r^-1 (sR -  eG)
		//               Q = r^-1 (sR + -eG)
		return this.g.mulAdd(s1, r, s2)
	}

	getKeyRecoveryParam(e, signature, Q, enc) {
		signature = Signature(signature, enc)
		if (signature.recoveryParam !== null) return signature.recoveryParam

		for (let i = 0; i < 4; i++) {
			let Qprime
			try {
				Qprime = this.recoverPubKey(e, signature, i, undefined)
			} catch (e) {
				continue
			}

			if (Qprime.eq(Q)) return i
		}
		throw new Error('Unable to find valid recovery factor')
	}
}
