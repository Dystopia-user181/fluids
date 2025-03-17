import { Vector3 } from "three";

export const worldBottom = -10;
const worldWidth = 6;
const restitution = 0;
const influence = 0.5;
function getSegmentKey(pos: Vector3) {
	return ((Math.round((pos.x + 30) / influence) << 8) +
		(Math.round((pos.y + 30) / influence) << 4) +
		Math.round((pos.z + 30) / influence)) & 4095;
}

const targetDensity = 0.7;

class FluidHandler {
	n = 0;
	r: Vector3[] = [];
	v: Vector3[] = [];
	segmentedR = Array(4096).fill(0).map(() => [] as number[]);

	densities: number[] = [];
	gradients: Vector3[] = [];
	pressureForceStrength = 0.005;

	constructor(number: number) {
		this.n = number;
		for (let i = 0; i < this.n; i++) {
			this.r.push(new Vector3(Math.random() * 8 - 4, Math.random() * 8 - 4, 0));
			this.v.push(new Vector3(0, 0, 0));
			this.densities.push(0.5);
			this.gradients.push(new Vector3(0, 0, 0));
			this.segmentedR[getSegmentKey(this.r[i])].push(i);
		}
	}

	getAllInSegment(pos: Vector3) {
		return this.segmentedR[getSegmentKey(pos)];
	}

	gravity(i: number) {
		const r = this.r[i].clone();
		return r.multiplyScalar(-Math.min(0.03 / (r.length() ** 3), 1 / r.length()));
	}

	calcDensity(u: number, dt: number) {
		this.densities[u] = 0.5;
		const effectivePos = this.r[u].clone().add(this.v[u].clone().multiplyScalar(dt));
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					const newPos = new Vector3(
						effectivePos.x + i * influence,
						effectivePos.y + j * influence,
						effectivePos.z + k * influence
					);
					let v;
					for (v of this.segmentedR[getSegmentKey(newPos)]) {
						// eslint-disable-next-line max-depth
						if (v === u) continue;
						const dr = new Vector3(
							this.r[v].x - effectivePos.x,
							this.r[v].y - effectivePos.y,
							this.r[v].z - effectivePos.z,
						);
						const magDrSq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
						// eslint-disable-next-line max-depth
						if (magDrSq > 0.25) continue;
						const density = (0.5 - 2 * Math.sqrt(magDrSq) + 2 * magDrSq);
						this.densities[u] += density;
					}
				}
			}
		}
	}

	calcPressureGrad(u: number, dt: number) {
		const effectivePos = this.r[u].clone().add(this.v[u].clone().multiplyScalar(dt));
		// cole equation: p is proportional to (rho / rho_0)^7 - 1
		const u1 = this.densities[u] / targetDensity - 1;
		const u2 = u1 * u1;
		const densityMul = u2 * u2 * u2;
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					const newPos = new Vector3(
						effectivePos.x + i * influence,
						effectivePos.y + j * influence,
						effectivePos.z + k * influence
					);
					let v;
					for (v of this.segmentedR[getSegmentKey(newPos)]) {
						// eslint-disable-next-line max-depth
						if (v === u) continue;
						const dr = new Vector3(
							this.r[v].x - effectivePos.x,
							this.r[v].y - effectivePos.y,
							this.r[v].z - effectivePos.z,
						);
						const magDrSq = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
						// eslint-disable-next-line max-depth
						if (magDrSq > 0.25) continue;
						// This is equivalent to dr.normalise().mul(grad density)
						dr.multiplyScalar((4 - 2 / Math.sqrt(magDrSq)) * densityMul);
						this.gradients[u].add(dr);
						this.gradients[v].sub(dr);
					}
				}
			}
		}
	}

	tick(dt: number) {
		const _dt = Math.min(dt, 0.2);
		const n = Math.ceil(_dt / 0.01);

		for (let i = 0; i < n; i++) {
			// let prev = calledTimes;
			this._tick(_dt / n);
			// console.log(calledTimes - prev)
		}
	}

	_tick(dt: number) {
		this.gradients.map(() => new Vector3(0, 0, 0));
		for (let i = 0; i < this.n; i++) {
			this.calcDensity(i, dt);
		}
		for (let i = 0; i < this.n; i++) {
			this.calcPressureGrad(i, dt);
		}
		const testThing = new Vector3(0, 0, 0);
		for (const g of this.gradients) testThing.add(g);
		// console.log(testThing)
		for (let i = 0; i < this.n; i++) {
			this.v[i].add(this.gravity(i).multiplyScalar(dt));
			this.v[i].add(this.gradients[i].multiplyScalar(dt * this.pressureForceStrength));
			const prevR = getSegmentKey(this.r[i]);
			this.r[i].add(this.v[i]);
			if (this.r[i].y < worldBottom) {
				this.r[i].y = worldBottom;
				this.v[i].y = -this.v[i].y * restitution;
			}
			if (Math.abs(this.r[i].x) > worldWidth) {
				this.r[i].x = Math.sign(this.r[i].x) * worldWidth;
				this.v[i].x *= -0.4;
				/* this.v[i].x = this.v[i].x * (0.04 ** dt) -
					Math.sign(this.r[i].x) * (Math.abs(this.r[i].x) - worldWidth) * dt; */
			}
			if (Math.abs(this.r[i].y) > worldWidth) {
				this.r[i].y = Math.sign(this.r[i].y) * worldWidth;
				this.v[i].y *= -0.4;
			}
			if (Math.abs(this.r[i].z) > worldWidth) {
				this.r[i].z = Math.sign(this.r[i].z) * worldWidth;
				this.v[i].z *= -0.4;
			}
			const newR = getSegmentKey(this.r[i]);
			if (prevR !== newR) {
				this.segmentedR[prevR].splice(this.segmentedR[prevR].indexOf(i), 1);
				this.segmentedR[newR].push(i);
			}
			this.v[i].multiplyScalar(0.7 ** dt);
		}
	}
}

export const Simulation = new FluidHandler(1000);

window.Simulation = Simulation;