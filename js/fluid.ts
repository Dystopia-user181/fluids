import { Vector3 } from "three";

export const worldBottom = -10;
const worldWidth = 4;
const restitution = 0;
const influence = 0.5;
function getSegmentKey(pos: Vector3) {
	return ((Math.floor((pos.x + worldWidth + 30) / influence) << 14) +
		(Math.floor((pos.y + worldWidth + 30) / influence) << 7) +
		Math.floor((pos.z + worldWidth + 30) / influence)) & 2097151;
}

const targetDensity = 0.7;

class FluidHandler {
	n = 0;
	r: Vector3[] = [];
	v: Vector3[] = [];
	segmentedR = Array(2097152).fill(0).map(() => new Set<number>());

	densities: number[] = [];
	gradients: Vector3[] = [];
	pressureForceStrength = 0.01;

	constructor(number: number) {
		this.n = number;
		for (let i = 0; i < this.n; i++) {
			this.r.push(new Vector3(Math.random() * 4 - 2, Math.random() * 4 - 2, 0));
			this.v.push(new Vector3(0, 0, 0));
			this.densities.push(0.5);
			this.gradients.push(new Vector3(0, 0, 0));
			this.segmentedR[getSegmentKey(this.r[i])].add(i);
		}
	}

	getAllInSegment(pos: Vector3) {
		return this.segmentedR[getSegmentKey(pos)];
	}

	gravity(i: number) {
		return new Vector3(Math.sin(Date.now() / 1500) * 0.003, -0.01, 0);
		const r = this.r[i].clone();
		return r.multiplyScalar(-Math.min(0.1 / (r.length() ** 3), 10 / r.length()));
	}

	calcDensity(u: number, dt: number) {
		this.densities[u] = 0.5;
		const effectivePos = this.r[u].clone().add(this.v[u].clone().multiplyScalar(dt));
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					const newPos = effectivePos.clone().add(new Vector3(i, j, k).multiplyScalar(influence));
					for (const j of this.segmentedR[getSegmentKey(newPos)]) {
						const dr = this.r[j].clone().sub(effectivePos);
						const magDr = dr.length();
						// eslint-disable-next-line max-depth
						if (magDr <= 0.0001 || magDr > 0.5) continue;
						const density = (0.5 - 2 * magDr + 2 * magDr * magDr);
						this.densities[u] += density;
					}
				}
			}
		}
	}

	calcPressureGrad(u: number, dt: number) {
		const effectivePos = this.r[u].clone().add(this.v[u].clone().multiplyScalar(dt));
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					const newPos = effectivePos.clone().add(new Vector3(i, j, k).multiplyScalar(influence));
					for (const j of this.segmentedR[getSegmentKey(newPos)]) {
						const dr = this.r[j].clone().sub(effectivePos);
						const magDr = dr.length();
						// eslint-disable-next-line max-depth
						if (magDr === 0 || magDr > 0.5) continue;
						// cole equation: p is proportional to (rho / rho_0)^7 - 1
						dr.normalize().multiplyScalar(-2 + 4 * magDr);
						this.gradients[u].add(dr.multiplyScalar(
							Math.min(this.densities[u] / targetDensity - 1, 3) ** 3
						));
						this.gradients[j].add(dr.multiplyScalar(-1));
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
				this.segmentedR[prevR].delete(i);
				this.segmentedR[newR].add(i);
			}
			this.v[i].multiplyScalar(0.7 ** dt);
		}
	}
}

export const Simulation = new FluidHandler(500);

window.Simulation = Simulation;