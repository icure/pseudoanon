import { BN } from 'bn.js'
import { utils } from '../utils'
import { BaseCurve } from './base'
import { inherits } from 'node:util'

export const mont = MontCurve

function MontCurve(conf) {
	BaseCurve.call(this, 'mont', conf)

	this.a = new BN(conf.a, 16).toRed(this.red)
	this.b = new BN(conf.b, 16).toRed(this.red)
	this.i4 = new BN(4).toRed(this.red).redInvm()
	this.two = new BN(2).toRed(this.red)
	this.a24 = this.i4.redMul(this.a.redAdd(this.two))
}
inherits(MontCurve, BaseCurve)

MontCurve.prototype.validate = function validate(point) {
	const x = point.normalize().x
	const x2 = x.redSqr()
	const rhs = x2.redMul(x).redAdd(x2.redMul(this.a)).redAdd(x)
	const y = rhs.redSqrt()

	return y.redSqr().cmp(rhs) === 0
}

function Point(curve, x, z) {
	BaseCurve.BasePoint.call(this, curve, 'projective')
	if (x === null && z === null) {
		this.x = this.curve.one
		this.z = this.curve.zero
	} else {
		this.x = new BN(x, 16)
		this.z = new BN(z, 16)
		if (!this.x.red) this.x = this.x.toRed(this.curve.red)
		if (!this.z.red) this.z = this.z.toRed(this.curve.red)
	}
}
inherits(Point, BaseCurve.BasePoint)

MontCurve.prototype.decodePoint = function decodePoint(bytes, enc) {
	return this.point(utils.toArray(bytes, enc), 1)
}

MontCurve.prototype.point = function point(x, z) {
	return new Point(this, x, z)
}

MontCurve.prototype.pointFromJSON = function pointFromJSON(obj) {
	return Point.fromJSON(this, obj)
}

Point.prototype.precompute = function precompute() {
	// No-op
}

Point.prototype._encode = function _encode() {
	return this.getX().toArray('be', this.curve.p.byteLength())
}

Point.fromJSON = function fromJSON(curve, obj) {
	return new Point(curve, obj[0], obj[1] || curve.one)
}

Point.prototype.inspect = function inspect() {
	if (this.isInfinity()) return '<EC Point Infinity>'
	return '<EC Point x: ' + this.x.fromRed().toString(16, 2) + ' z: ' + this.z.fromRed().toString(16, 2) + '>'
}

Point.prototype.isInfinity = function isInfinity() {
	// XXX This code assumes that zero is always zero in red
	return this.z.cmpn(0) === 0
}

Point.prototype.dbl = function dbl() {
	// http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#doubling-dbl-1987-m-3
	// 2M + 2S + 4A

	// A = X1 + Z1
	const a = this.x.redAdd(this.z)
	// AA = A^2
	const aa = a.redSqr()
	// B = X1 - Z1
	const b = this.x.redSub(this.z)
	// BB = B^2
	const bb = b.redSqr()
	// C = AA - BB
	const c = aa.redSub(bb)
	// X3 = AA * BB
	const nx = aa.redMul(bb)
	// Z3 = C * (BB + A24 * C)
	const nz = c.redMul(bb.redAdd(this.curve.a24.redMul(c)))
	return this.curve.point(nx, nz)
}

Point.prototype.add = function add() {
	throw new Error('Not supported on Montgomery curve')
}

Point.prototype.diffAdd = function diffAdd(p, diff) {
	// http://hyperelliptic.org/EFD/g1p/auto-montgom-xz.html#diffadd-dadd-1987-m-3
	// 4M + 2S + 6A

	// A = X2 + Z2
	const a = this.x.redAdd(this.z)
	// B = X2 - Z2
	const b = this.x.redSub(this.z)
	// C = X3 + Z3
	const c = p.x.redAdd(p.z)
	// D = X3 - Z3
	const d = p.x.redSub(p.z)
	// DA = D * A
	const da = d.redMul(a)
	// CB = C * B
	const cb = c.redMul(b)
	// X5 = Z1 * (DA + CB)^2
	const nx = diff.z.redMul(da.redAdd(cb).redSqr())
	// Z5 = X1 * (DA - CB)^2
	const nz = diff.x.redMul(da.redISub(cb).redSqr())
	return this.curve.point(nx, nz)
}

Point.prototype.mul = function mul(k) {
	const t = k.clone()
	// eslint-disable-next-line @typescript-eslint/no-this-alias
	let a = this // (N / 2) * Q + Q
	let b = this.curve.point(null, null) // (N / 2) * Q
	// eslint-disable-next-line @typescript-eslint/no-this-alias
	const c = this // Q

	let bits

	for (bits = []; t.cmpn(0) !== 0; t.iushrn(1)) bits.push(t.andln(1))

	for (let i = bits.length - 1; i >= 0; i--) {
		if (bits[i] === 0) {
			// N * Q + Q = ((N / 2) * Q + Q)) + (N / 2) * Q
			a = a.diffAdd(b, c)
			// N * Q = 2 * ((N / 2) * Q + Q))
			b = b.dbl()
		} else {
			// N * Q = ((N / 2) * Q + Q) + ((N / 2) * Q)
			b = a.diffAdd(b, c)
			// N * Q + Q = 2 * ((N / 2) * Q + Q)
			a = a.dbl()
		}
	}
	return b
}

Point.prototype.mulAdd = function mulAdd() {
	throw new Error('Not supported on Montgomery curve')
}

Point.prototype.jumlAdd = function jumlAdd() {
	throw new Error('Not supported on Montgomery curve')
}

Point.prototype.eq = function eq(other) {
	return this.getX().cmp(other.getX()) === 0
}

Point.prototype.normalize = function normalize() {
	this.x = this.x.redMul(this.z.redInvm())
	this.z = this.curve.one
	return this
}

Point.prototype.getX = function getX() {
	// Normalize coordinates
	this.normalize()

	return this.x.fromRed()
}
