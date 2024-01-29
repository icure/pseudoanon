import { EC } from './lib/elliptic/ec'

const samples = [
	[
		'MQ==',
		8,
		'MQEAAAAAAAAAAA==',
		'wAIaELoSge8ZH9YhlIeahODkJDe23j1NUsQqwm32j1o+CdG6lIKhi1KqdONLRsXh+ciZjddZtm2dShPpp3K5aDxD8',
	],
	[
		'MTI=',
		8,
		'MTICAAAAAAAAAAA=',
		'AURV0WACmddkLIQK2IF76J0XGOygU+mSZ5gnbQFQ7WSyEUq9H/c0738e2pwVTiSxcI1xz1dGqTcYfotir3K0LQrz',
	],
	[
		'MTIz',
		8,
		'MTIzAwAAAAAAAAAC',
		'QRJBFXWQueauTw1qeNzTM+3qZtLLbykEN8LCmrn29Nv3IXHlLFa9Hhd/iuRBHDsgv4usav38XCWgxroQ8zWruQ4=',
	],
] as [string, number, string, string][]

const curveP521 = new EC('p521')

const identifierToECPoint = (identifier: string, bufferSize: number) => {
	const idBytes = Buffer.from(identifier, 'base64')
	const xBufferSize = bufferSize
	const xBytes = Buffer.alloc(1 + idBytes.length + 1 + xBufferSize, 0)
	//xBytes[0] = xBytes[0] + 1; // added in draft cookbook 1.4, but does not give correct results

	idBytes.copy(xBytes, 1)
	xBytes[1 + idBytes.length] = idBytes.length //added in draft cookbook 1.4

	let point

	//console.log(this.curveP521.a, this.curveP521.b, this.curveP521.curve.p.toString(16));

	const xCoordinate = curveP521.keyFromPrivate(xBytes, undefined).getPrivate(undefined) // - same result as below
	//xCoordinate = new BN(xBytes);
	const maxiterations = 32000
	let i = 0
	do {
		i++
		try {
			point = curveP521.curve.pointFromX(xCoordinate, true)
			if (point) {
				point.getY()
				break
			}
		} catch (error) {
			point = null
			xCoordinate.iaddn(1)
		}
	} while (!point && i < maxiterations)

	return point
}

samples.forEach(([identifier, bufferSize, xc, yc]) => {
	const idBytes = Buffer.from(identifier, 'base64')
	const xBufferSize = bufferSize
	const xBytes = Buffer.alloc(1 + idBytes.length + 1 + xBufferSize, 0)
	idBytes.copy(xBytes, 1)
	xBytes[0] = 1 //added in draft cookbook 1.4
	xBytes[1 + idBytes.length] = idBytes.length //added in draft cookbook 1.4

	samples.forEach(([identifier, bufferSize, xc, yc]) => {
		const ecPoint = identifierToECPoint(identifier, bufferSize)

		const x = Buffer.from(ecPoint.getX().toArray()).toString('base64')
		const y = Buffer.from(ecPoint.getY().toArray()).toString('base64')

		const okx = xc === x
		const oky = yc === y

		console.log('identifier', identifier, '\nx', x, '\nokx', okx, '\ny', y, '\noky', oky)
		return { identifier: identifier, ecPoint: ecPoint, x: x, y: y, oky: oky }
	})
})
