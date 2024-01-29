import { BN } from 'bn.js'
import { utils } from '../utils'
import { BaseCurve } from './base'
import { inherits } from 'node:util'

const assert = utils.assert
export const edwards = EdwardsCurve

function EdwardsCurve(conf) {
	// NOTE: Important as we are creating point in BaseCurve.call()
	this.twisted = (conf.a | 0) !== 1
	this.mOneA = this.twisted && (conf.a | 0) === -1
	this.extended = this.mOneA

	BaseCurve.call(this, 'edwards', conf)

	this.a = new BN(conf.a, 16).umod(this.red.m)
	this.a = this.a.toRed(this.red)
	this.c = new BN(conf.c, 16).toRed(this.red)
	this.c2 = this.c.redSqr()
	this.d = new BN(conf.d, 16).toRed(this.red)
	this.dd = this.d.redAdd(this.d)

	assert(!this.twisted || this.c.fromRed().cmpn(1) === 0)
	this.oneC = (conf.c | 0) === 1
}
inherits(EdwardsCurve, BaseCurve)

EdwardsCurve.prototype._mulA = function _mulA(num) {
	if (this.mOneA) return num.redNeg()
	else return this.a.redMul(num)
}

EdwardsCurve.prototype._mulC = function _mulC(num) {
	if (this.oneC) return num
	else return this.c.redMul(num)
}

// Just for compatibility with Short curve
EdwardsCurve.prototype.jpoint = function jpoint(x, y, z, t) {
	return this.point(x, y, z, t)
}

EdwardsCurve.prototype.pointFromX = function pointFromX(x, odd) {
	x = new BN(x, 16)
	if (!x.red) x = x.toRed(this.red)

	const x2 = x.redSqr()
	const rhs = this.c2.redSub(this.a.redMul(x2))
	const lhs = this.one.redSub(this.c2.redMul(this.d).redMul(x2))

	const y2 = rhs.redMul(lhs.redInvm())
	let y = y2.redSqrt()
	if (y.redSqr().redSub(y2).cmp(this.zero) !== 0) throw new Error('invalid point')

	const isOdd = y.fromRed().isOdd()
	if ((odd && !isOdd) || (!odd && isOdd)) y = y.redNeg()

	return this.point(x, y)
}

EdwardsCurve.prototype.pointFromY = function pointFromY(y, odd) {
	y = new BN(y, 16)
	if (!y.red) y = y.toRed(this.red)

	// x^2 = (y^2 - c^2) / (c^2 d y^2 - a)
	const y2 = y.redSqr()
	const lhs = y2.redSub(this.c2)
	const rhs = y2.redMul(this.d).redMul(this.c2).redSub(this.a)
	const x2 = lhs.redMul(rhs.redInvm())

	if (x2.cmp(this.zero) === 0) {
		if (odd) throw new Error('invalid point')
		else return this.point(this.zero, y)
	}

	let x = x2.redSqrt()
	if (x.redSqr().redSub(x2).cmp(this.zero) !== 0) throw new Error('invalid point')

	if (x.fromRed().isOdd() !== odd) x = x.redNeg()

	return this.point(x, y)
}

EdwardsCurve.prototype.validate = function validate(point) {
	if (point.isInfinity()) return true

	// Curve: A * X^2 + Y^2 = C^2 * (1 + D * X^2 * Y^2)
	point.normalize()

	const x2 = point.x.redSqr()
	const y2 = point.y.redSqr()
	const lhs = x2.redMul(this.a).redAdd(y2)
	const rhs = this.c2.redMul(this.one.redAdd(this.d.redMul(x2).redMul(y2)))

	return lhs.cmp(rhs) === 0
}

function Point(curve, x, y, z, t) {
	BaseCurve.BasePoint.call(this, curve, 'projective')
	if (x === null && y === null && z === null) {
		this.x = this.curve.zero
		this.y = this.curve.one
		this.z = this.curve.one
		this.t = this.curve.zero
		this.zOne = true
	} else {
		this.x = new BN(x, 16)
		this.y = new BN(y, 16)
		this.z = z ? new BN(z, 16) : this.curve.one
		this.t = t && new BN(t, 16)
		if (!this.x.red) this.x = this.x.toRed(this.curve.red)
		if (!this.y.red) this.y = this.y.toRed(this.curve.red)
		if (!this.z.red) this.z = this.z.toRed(this.curve.red)
		if (this.t && !this.t.red) this.t = this.t.toRed(this.curve.red)
		this.zOne = this.z === this.curve.one

		// Use extended coordinates
		if (this.curve.extended && !this.t) {
			this.t = this.x.redMul(this.y)
			if (!this.zOne) this.t = this.t.redMul(this.z.redInvm())
		}
	}
}
inherits(Point, BaseCurve.BasePoint)

EdwardsCurve.prototype.pointFromJSON = function pointFromJSON(obj) {
	return Point.fromJSON(this, obj)
}

EdwardsCurve.prototype.point = function point(x, y, z, t) {
	return new Point(this, x, y, z, t)
}

Point.fromJSON = function fromJSON(curve, obj) {
	return new Point(curve, obj[0], obj[1], obj[2], undefined)
}

Point.prototype.inspect = function inspect() {
	if (this.isInfinity()) return '<EC Point Infinity>'
	return (
		'<EC Point x: ' +
		this.x.fromRed().toString(16, 2) +
		' y: ' +
		this.y.fromRed().toString(16, 2) +
		' z: ' +
		this.z.fromRed().toString(16, 2) +
		'>'
	)
}

Point.prototype.isInfinity = function isInfinity() {
	// XXX This code assumes that zero is always zero in red
	return this.x.cmpn(0) === 0 && (this.y.cmp(this.z) === 0 || (this.zOne && this.y.cmp(this.curve.c) === 0))
}

Point.prototype._extDbl = function _extDbl() {
	// hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html
	//     #doubling-dbl-2008-hwcd
	// 4M + 4S

	// A = X1^2
	const a = this.x.redSqr()
	// B = Y1^2
	const b = this.y.redSqr()
	// C = 2 * Z1^2
	let c = this.z.redSqr()
	c = c.redIAdd(c)
	// D = a * A
	const d = this.curve._mulA(a)
	// E = (X1 + Y1)^2 - A - B
	const e = this.x.redAdd(this.y).redSqr().redISub(a).redISub(b)
	// G = D + B
	const g = d.redAdd(b)
	// F = G - C
	const f = g.redSub(c)
	// H = D - B
	const h = d.redSub(b)
	// X3 = E * F
	const nx = e.redMul(f)
	// Y3 = G * H
	const ny = g.redMul(h)
	// T3 = E * H
	const nt = e.redMul(h)
	// Z3 = F * G
	const nz = f.redMul(g)
	return this.curve.point(nx, ny, nz, nt)
}

Point.prototype._projDbl = function _projDbl() {
	// hyperelliptic.org/EFD/g1p/auto-twisted-projective.html
	//     #doubling-dbl-2008-bbjlp
	//     #doubling-dbl-2007-bl
	// and others
	// Generally 3M + 4S or 2M + 4S

	// B = (X1 + Y1)^2
	const b = this.x.redAdd(this.y).redSqr()
	// C = X1^2
	const c = this.x.redSqr()
	// D = Y1^2
	const d = this.y.redSqr()

	let nx
	let ny
	let nz
	let e
	let h
	let j
	if (this.curve.twisted) {
		// E = a * C
		e = this.curve._mulA(c)
		// F = E + D
		const f = e.redAdd(d)
		if (this.zOne) {
			// X3 = (B - C - D) * (F - 2)
			nx = b.redSub(c).redSub(d).redMul(f.redSub(this.curve.two))
			// Y3 = F * (E - D)
			ny = f.redMul(e.redSub(d))
			// Z3 = F^2 - 2 * F
			nz = f.redSqr().redSub(f).redSub(f)
		} else {
			// H = Z1^2
			h = this.z.redSqr()
			// J = F - 2 * H
			j = f.redSub(h).redISub(h)
			// X3 = (B-C-D)*J
			nx = b.redSub(c).redISub(d).redMul(j)
			// Y3 = F * (E - D)
			ny = f.redMul(e.redSub(d))
			// Z3 = F * J
			nz = f.redMul(j)
		}
	} else {
		// E = C + D
		e = c.redAdd(d)
		// H = (c * Z1)^2
		h = this.curve._mulC(this.z).redSqr()
		// J = E - 2 * H
		j = e.redSub(h).redSub(h)
		// X3 = c * (B - E) * J
		nx = this.curve._mulC(b.redISub(e)).redMul(j)
		// Y3 = c * E * (C - D)
		ny = this.curve._mulC(e).redMul(c.redISub(d))
		// Z3 = E * J
		nz = e.redMul(j)
	}
	return this.curve.point(nx, ny, nz)
}

Point.prototype.dbl = function dbl() {
	if (this.isInfinity()) return this

	// Double in extended coordinates
	if (this.curve.extended) return this._extDbl()
	else return this._projDbl()
}

Point.prototype._extAdd = function _extAdd(p) {
	// hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html
	//     #addition-add-2008-hwcd-3
	// 8M

	// A = (Y1 - X1) * (Y2 - X2)
	const a = this.y.redSub(this.x).redMul(p.y.redSub(p.x))
	// B = (Y1 + X1) * (Y2 + X2)
	const b = this.y.redAdd(this.x).redMul(p.y.redAdd(p.x))
	// C = T1 * k * T2
	const c = this.t.redMul(this.curve.dd).redMul(p.t)
	// D = Z1 * 2 * Z2
	const d = this.z.redMul(p.z.redAdd(p.z))
	// E = B - A
	const e = b.redSub(a)
	// F = D - C
	const f = d.redSub(c)
	// G = D + C
	const g = d.redAdd(c)
	// H = B + A
	const h = b.redAdd(a)
	// X3 = E * F
	const nx = e.redMul(f)
	// Y3 = G * H
	const ny = g.redMul(h)
	// T3 = E * H
	const nt = e.redMul(h)
	// Z3 = F * G
	const nz = f.redMul(g)
	return this.curve.point(nx, ny, nz, nt)
}

Point.prototype._projAdd = function _projAdd(p) {
	// hyperelliptic.org/EFD/g1p/auto-twisted-projective.html
	//     #addition-add-2008-bbjlp
	//     #addition-add-2007-bl
	// 10M + 1S

	// A = Z1 * Z2
	const a = this.z.redMul(p.z)
	// B = A^2
	const b = a.redSqr()
	// C = X1 * X2
	const c = this.x.redMul(p.x)
	// D = Y1 * Y2
	const d = this.y.redMul(p.y)
	// E = d * C * D
	const e = this.curve.d.redMul(c).redMul(d)
	// F = B - E
	const f = b.redSub(e)
	// G = B + E
	const g = b.redAdd(e)
	// X3 = A * F * ((X1 + Y1) * (X2 + Y2) - C - D)
	const tmp = this.x.redAdd(this.y).redMul(p.x.redAdd(p.y)).redISub(c).redISub(d)
	const nx = a.redMul(f).redMul(tmp)
	let ny
	let nz
	if (this.curve.twisted) {
		// Y3 = A * G * (D - a * C)
		ny = a.redMul(g).redMul(d.redSub(this.curve._mulA(c)))
		// Z3 = F * G
		nz = f.redMul(g)
	} else {
		// Y3 = A * G * (D - C)
		ny = a.redMul(g).redMul(d.redSub(c))
		// Z3 = c * F * G
		nz = this.curve._mulC(f).redMul(g)
	}
	return this.curve.point(nx, ny, nz)
}

Point.prototype.add = function add(p) {
	if (this.isInfinity()) return p
	if (p.isInfinity()) return this

	if (this.curve.extended) return this._extAdd(p)
	else return this._projAdd(p)
}

Point.prototype.mul = function mul(k) {
	if (this._hasDoubles(k)) return this.curve._fixedNafMul(this, k)
	else return this.curve._wnafMul(this, k)
}

Point.prototype.mulAdd = function mulAdd(k1, p, k2) {
	return this.curve._wnafMulAdd(1, [this, p], [k1, k2], 2, false)
}

Point.prototype.jmulAdd = function jmulAdd(k1, p, k2) {
	return this.curve._wnafMulAdd(1, [this, p], [k1, k2], 2, true)
}

Point.prototype.normalize = function normalize() {
	if (this.zOne) return this

	// Normalize coordinates
	const zi = this.z.redInvm()
	this.x = this.x.redMul(zi)
	this.y = this.y.redMul(zi)
	if (this.t) this.t = this.t.redMul(zi)
	this.z = this.curve.one
	this.zOne = true
	return this
}

Point.prototype.neg = function neg() {
	return this.curve.point(this.x.redNeg(), this.y, this.z, this.t && this.t.redNeg())
}

Point.prototype.getX = function getX() {
	this.normalize()
	return this.x.fromRed()
}

Point.prototype.getY = function getY() {
	this.normalize()
	return this.y.fromRed()
}

Point.prototype.eq = function eq(other) {
	return this === other || (this.getX().cmp(other.getX()) === 0 && this.getY().cmp(other.getY()) === 0)
}

Point.prototype.eqXToP = function eqXToP(x) {
	const rx = x.toRed(this.curve.red).redMul(this.z)
	if (this.x.cmp(rx) === 0) return true

	const xc = x.clone()
	const t = this.curve.redN.redMul(this.z)
	for (;;) {
		xc.iadd(this.curve.n)
		if (xc.cmp(this.curve.p) >= 0) return false

		rx.redIAdd(t)
		if (this.x.cmp(rx) === 0) return true
	}
}

// Compatibility with BaseCurve
Point.prototype.toP = Point.prototype.normalize
Point.prototype.mixedAdd = Point.prototype.add
