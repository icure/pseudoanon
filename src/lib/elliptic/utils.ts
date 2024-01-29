import { BN } from 'bn.js'
import { assert } from 'minimalistic-assert'
import { toArray, zero2, toHex, encode } from 'minimalistic-crypto-utils'

export const utils = {
	assert: (val: any, msg?: string) => {
		if (!val) throw new Error(msg || 'Assertion failed')
	},
	toArray: toArray,
	zero2: zero2,
	toHex: toHex,
	encode: encode,
	getNAF: getNAF,
	cachedProperty: cachedProperty,
	parseBytes: parseBytes,
	getJSF: getJSF,
	intFromLE: intFromLE,
}

// Represent num in a w-NAF form
function getNAF(num, w, bits) {
	const naf = new Array(Math.max(num.bitLength(), bits) + 1)
	naf.fill(0)

	const ws = 1 << (w + 1)
	const k = num.clone()

	for (let i = 0; i < naf.length; i++) {
		let z
		const mod = k.andln(ws - 1)
		if (k.isOdd()) {
			if (mod > (ws >> 1) - 1) z = (ws >> 1) - mod
			else z = mod
			k.isubn(z)
		} else {
			z = 0
		}

		naf[i] = z
		k.iushrn(1)
	}

	return naf
}

// Represent k1, k2 in a Joint Sparse Form
function getJSF(k1, k2) {
	const jsf = [[], []]

	k1 = k1.clone()
	k2 = k2.clone()
	let d1 = 0
	let d2 = 0
	let m8
	while (k1.cmpn(-d1) > 0 || k2.cmpn(-d2) > 0) {
		// First phase
		let m14 = (k1.andln(3) + d1) & 3
		let m24 = (k2.andln(3) + d2) & 3
		if (m14 === 3) m14 = -1
		if (m24 === 3) m24 = -1
		let u1
		if ((m14 & 1) === 0) {
			u1 = 0
		} else {
			m8 = (k1.andln(7) + d1) & 7
			if ((m8 === 3 || m8 === 5) && m24 === 2) u1 = -m14
			else u1 = m14
		}
		jsf[0].push(u1)

		let u2
		if ((m24 & 1) === 0) {
			u2 = 0
		} else {
			m8 = (k2.andln(7) + d2) & 7
			if ((m8 === 3 || m8 === 5) && m14 === 2) u2 = -m24
			else u2 = m24
		}
		jsf[1].push(u2)

		// Second phase
		if (2 * d1 === u1 + 1) d1 = 1 - d1
		if (2 * d2 === u2 + 1) d2 = 1 - d2
		k1.iushrn(1)
		k2.iushrn(1)
	}

	return jsf
}

function cachedProperty(obj, name, computer) {
	const key = '_' + name
	obj.prototype[name] = function cachedProperty() {
		return this[key] !== undefined ? this[key] : (this[key] = computer.call(this))
	}
}

function parseBytes(bytes) {
	return typeof bytes === 'string' ? utils.toArray(bytes, 'hex') : bytes
}

function intFromLE(bytes) {
	return new BN(bytes, 'hex', 'le')
}
