import { Vector3 } from "three";

export const worldBottom = -10;
const worldWidth = 6;
const restitution = 0.7;
const influence = 0.5;
function getSegmentKey(pos: Vector3) {
	return ((Math.round((pos.x + 30) / influence) << 8) +
		(Math.round((pos.y + 30) / influence) << 4) +
		Math.round((pos.z + 30) / influence)) & 4095;
}

const targetDensity = 0.75, selfDensity = 0.5;

class FluidHandler {
	n = 0;
	r: Vector3[] = [];
	v: Vector3[] = [];
	segmentedR = Array(4096).fill(0).map(() => [] as number[]);

	densities: number[] = [];
	gradP: Vector3[] = [];
	lapV: Vector3[] = [];
	forces = {
		mag: [] as Vector3[],
		by: [] as number[],
		to: [] as number[],
	};

	pressureForceStrength = 3;
	viscosity = 10;

	timeElapsed = 0;
	interactionCounter = 0;

	constructor(number: number) {
		this.n = number;
		for (let i = 0; i < this.n; i++) {
			this.r.push(new Vector3(Math.random() * 12 - 6, Math.random() * 12 - 6, 0));
			this.v.push(new Vector3(0, 0, 0));
			this.densities.push(selfDensity);
			this.gradP.push(new Vector3(0, 0, 0));
			this.lapV.push(new Vector3(0, 0, 0));
			this.segmentedR[getSegmentKey(this.r[i])].push(i);
		}
	}

	getAllInSegment(pos: Vector3) {
		return this.segmentedR[getSegmentKey(pos)];
	}

	gravity(i: number) {
		return new Vector3(0, -2, 0);
		const lengthSq = this.r[i].x * this.r[i].x + (this.r[i].y + worldWidth) * (this.r[i].y + worldWidth);
		const convection = new Vector3(0, 20 * Math.exp(-lengthSq), 0);
		return new Vector3(0, -2, 0).add(convection);
		const phase = this.timeElapsed * 0.0609;
		const r = this.r[i].clone().add(new Vector3(2 * Math.cos(phase), 2 * Math.sin(phase), 0));
		const g1 = r.clone().multiplyScalar(-Math.min(0.1 / (r.lengthSq() ** 1.5), 1));
		const r2 = r.sub(new Vector3(4 * Math.cos(phase), 4 * Math.sin(phase), 0));
		const g2 = r2.clone().multiplyScalar(-Math.min(0.1 / (r2.lengthSq() ** 1.5), 1));
		return g1.add(g2);
	}

	calcDensity(u: number) {
		this.densities[u] = selfDensity;
		const start = this.forces.mag.length;
		const lapV = this.lapV[u];
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					const newPos = new Vector3(
						this.r[u].x + i * influence,
						this.r[u].y + j * influence,
						this.r[u].z + k * influence
					);
					let v;
					for (v of this.segmentedR[getSegmentKey(newPos)]) {
						// eslint-disable-next-line max-depth
						if (v === u) continue;
						const dr = new Vector3(
							this.r[v].x - this.r[u].x,
							this.r[v].y - this.r[u].y,
							this.r[v].z - this.r[u].z,
						);
						const magDrSq = dr.lengthSq();
						// eslint-disable-next-line max-depth
						if (magDrSq < 1e-15 || magDrSq > 0.25) continue;
						this.interactionCounter++;
						const density = (0.5 - 2 * Math.sqrt(magDrSq) + 2 * magDrSq);
						this.densities[u] += density;
						this.forces.mag.push(dr.multiplyScalar(4 - 2 / Math.sqrt(magDrSq)));
						this.forces.to.push(u);
						this.forces.by.push(v);
						const dv = new Vector3(
							this.v[v].x - this.v[u].x,
							this.v[v].y - this.v[u].y,
							this.v[v].z - this.v[u].z,
						);
						lapV.add(dv.multiplyScalar(density));
					}
				}
			}
		}
		const u1 = (this.densities[u] - targetDensity) / targetDensity;
		const u2 = u1 * u1;
		const densityMultiplier = Math.min(u2 * u2 * u2, 10) * this.pressureForceStrength;
		const end = this.forces.mag.length;
		for (let i = start; i < end; i++) {
			this.forces.mag[i].multiplyScalar(densityMultiplier);
		}
	}

	finalisePressureGrads() {
		const nForces = this.forces.mag.length;
		for (let i = 0; i < nForces; i++) {
			this.gradP[this.forces.to[i]].add(this.forces.mag[i]);
			this.gradP[this.forces.by[i]].sub(this.forces.mag[i]);
		}
		this.forces.mag = [];
		this.forces.by = [];
		this.forces.to = [];
	}

	tick(dt: number) {
		const _dt = Math.min(dt, 0.02);
		const n = Math.ceil(_dt / 1);

		for (let i = 0; i < n; i++) {
			// let prev = calledTimes;
			this._tick(_dt / n);
			// console.log(calledTimes - prev)
		}
	}

	_tick(dt: number) {
		this.interactionCounter = 0;
		this.timeElapsed += dt;
		this.gradP = this.gradP.map(() => new Vector3(0, 0, 0));
		this.lapV = this.lapV.map(() => new Vector3(0, 0, 0));
		for (let i = 0; i < this.n; i++) {
			this.calcDensity(i);
		}
		this.finalisePressureGrads();
		for (let i = 0; i < this.n; i++) {
			this.v[i].add(this.gravity(i).multiplyScalar(dt));
			this.v[i].add(this.gradP[i].clone().multiplyScalar(dt));
			this.v[i].add(this.lapV[i].clone().multiplyScalar(dt * this.viscosity));
			const prevR = getSegmentKey(this.r[i]);
			this.r[i].add(this.v[i].clone().multiplyScalar(dt));
			if (this.r[i].y < worldBottom) {
				this.r[i].y = worldBottom;
				this.v[i].y = -this.v[i].y * restitution;
			}
			if (Math.abs(this.r[i].x) > worldWidth) {
				this.r[i].x = Math.sign(this.r[i].x) * 2 * worldWidth - this.r[i].x;
				this.v[i].x *= -restitution;
				/* this.v[i].x = this.v[i].x * (0.04 ** dt) -
					Math.sign(this.r[i].x) * (Math.abs(this.r[i].x) - worldWidth) * dt; */
			}
			if (Math.abs(this.r[i].y) > worldWidth) {
				this.r[i].y = Math.sign(this.r[i].y) * 2 * worldWidth - this.r[i].y;
				this.v[i].y *= -restitution;
			}
			if (Math.abs(this.r[i].z) > worldWidth) {
				this.r[i].z = Math.sign(this.r[i].z) * 2 * worldWidth - this.r[i].z;
				this.v[i].z *= -restitution;
			}
			const newR = getSegmentKey(this.r[i]);
			if (prevR !== newR) {
				this.segmentedR[prevR].splice(this.segmentedR[prevR].indexOf(i), 1);
				this.segmentedR[newR].push(i);
			}
		}
	}
}

export const Simulation = new FluidHandler(3000);

// @ts-ignore
window.Simulation = Simulation;