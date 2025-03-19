import { Vector3 } from "three";

export const worldBottom = -10;
const worldWidth = 7;
const restitution = 0.1;
const influence = 0.5;
function interLeave3(input: number) {
	let x = input;
	x = (x | x << 8) & 0xf00f00f00f00f;
	x = (x | x << 4) & 0x30c30c30c30c3;
	x = (x | x << 2) & 0x0249249249249;
	return x;
}
function getSegmentKey(pos: Vector3) {
	return (interLeave3(Math.round((pos.x + 10) / influence)) +
		(interLeave3(Math.round((pos.y + 10) / influence)) << 1) +
		(interLeave3(Math.round((pos.z + 10) / influence)) << 2)) & 4095;
}

const targetDensity = 0.75, selfDensity = 0.5;

class FluidHandler {
	n = 0;
	r: Vector3[] = [];
	v: Vector3[] = [];
	chunkedR = Array(4096).fill(0).map(() => new Set<number>());
	inChunk: number[] = [];

	densities: number[] = [];
	gradP: Vector3[] = [];
	lapV: Vector3[] = [];

	pressureForceStrength = 3;
	viscosity = 2;
	tension = 3;

	timeElapsed = 0;
	interactionCounter = 0;

	constructor(number: number) {
		this.n = number;
		for (let i = 0; i < this.n; i++) {
			this.r.push(new Vector3(Math.random() * 7 - 7, Math.random() * 7 - 7, 0));
			this.v.push(new Vector3(0, 0, 0));
			this.densities.push(selfDensity);
			this.gradP.push(new Vector3(0, 0, 0));
			this.lapV.push(new Vector3(0, 0, 0));
			this.inChunk.push(4098);
			this.updateChunk(i);
		}
	}

	updateChunk(u: number, _prevPos?: Vector3) {
		const newKey = getSegmentKey(this.r[u]);
		if (this.inChunk[u] === newKey && !_prevPos) return;
		const prevPos = _prevPos ? _prevPos : new Vector3(Infinity, Infinity, Infinity);
		const dx = Math.round(this.r[u].x / influence) - Math.round(prevPos.x / influence);
		const dy = Math.round(this.r[u].y / influence) - Math.round(prevPos.y / influence);
		const dz = Math.round(this.r[u].z / influence) - Math.round(prevPos.z / influence);
		if (_prevPos) {
			for (let i = -1; i < 2; i++) {
				for (let j = -1; j < 2; j++) {
					for (let k = -1; k < 2; k++) {
						// eslint-disable-next-line max-depth
						if (-dx + i >= -1 && -dx + i <= 1 &&
							-dy + j >= -1 && -dy + j <= 1 &&
							-dz + k >= -1 && -dz + k <= 1) continue;
						const newPos = new Vector3(
							prevPos.x + i * influence,
							prevPos.y + j * influence,
							prevPos.z + k * influence
						);
						const key = getSegmentKey(newPos);
						this.chunkedR[key].delete(u);
					}
				}
			}
		}
		this.inChunk[u] = newKey;
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					// eslint-disable-next-line max-depth
					if (dx + i >= -1 && dx + i <= 1 &&
						dy + j >= -1 && dy + j <= 1 &&
						dz + k >= -1 && dz + k <= 1) continue;
					const newPos = new Vector3(
						this.r[u].x + i * influence,
						this.r[u].y + j * influence,
						this.r[u].z + k * influence
					);
					const key = getSegmentKey(newPos);
					this.chunkedR[key].add(u);
				}
			}
		}
	}

	gravity(i: number) {
		//return new Vector3(0, -2, 0);
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

	calcProperties() {
		const sortedIds = Array.from({ length: this.n }, (_, i) => i).sort((a, b) => this.inChunk[a] - this.inChunk[b]);
		const corrChunk = sortedIds.map(x => this.inChunk[x]);
		let ptr = 0;
		while (ptr < this.n) {
			const chunkBegins = ptr;
			const chunk = corrChunk[ptr];
			while (corrChunk[ptr] === corrChunk[ptr + 1]) ptr++;
			ptr++;
			const chunkEnds = ptr;
			// This is to help with caching
			const corrR = sortedIds.slice(chunkBegins, chunkEnds).map(x => this.r[x].clone());
			const corrV = sortedIds.slice(chunkBegins, chunkEnds).map(x => this.v[x].clone());
			const lapV = corrR.map(() => new Vector3());
			const forcesMag = [];
			const forcesTo = [];
			const forcesBy = [];
			for (const v of this.chunkedR[chunk]) {
				const { x: r1x, y: r1y, z: r1z } = this.r[v];
				const { x: v1x, y: v1y, z: v1z } = this.v[v];
				for (let i = chunkBegins; i < chunkEnds; i++) {
					const u = sortedIds[i];
					if (v === u) continue;
					// Putting const dr = ... After the calculation seemed inconsequential but I was able to
					// improve performance slightly with it from ~54fps to 60fps (when not running the renderer)
					// so I'm not questioning it
					const { x: r2x, y: r2y, z: r2z } = corrR[i - chunkBegins];
					const drx = r1x - r2x;
					const dry = r1y - r2y;
					const drz = r1z - r2z;
					const magDrSq = drx * drx + dry * dry + drz * drz;
					if (magDrSq < 1e-15 || magDrSq > 0.25) continue;
					const dr = new Vector3(
						drx,
						dry,
						drz,
					);
					this.gradP[u].add(dr.clone().multiplyScalar(this.tension));
					this.gradP[v].sub(dr.clone().multiplyScalar(this.tension));
					this.interactionCounter++;
					const density = (0.5 - 2 * Math.sqrt(magDrSq) + 2 * magDrSq);
					this.densities[u] += density;
					forcesMag.push(dr.multiplyScalar(4 - 2 / Math.sqrt(magDrSq)));
					forcesTo.push(u);
					forcesBy.push(v);
					const { x: v2x, y: v2y, z: v2z } = corrV[i - chunkBegins];
					const dv = new Vector3(
						v1x - v2x,
						v1y - v2y,
						v1z - v2z,
					);
					lapV[i - chunkBegins].add(dv.multiplyScalar(density));
				}
			}
			for (let i = chunkBegins; i < chunkEnds; i++) {
				const u = sortedIds[i];
				this.lapV[u].copy(lapV[i - chunkBegins]);
			}
			for (let i = 0; i < forcesMag.length; i++) {
				const u = forcesTo[i];
				const u1 = (this.densities[u]) / targetDensity;
				const u2 = u1 * u1;
				const densityMultiplier = Math.min(u2 * u2 * u2, 10) * this.pressureForceStrength;
				forcesMag[i].multiplyScalar(densityMultiplier);
				this.gradP[u].add(forcesMag[i]);
				this.gradP[forcesBy[i]].sub(forcesMag[i]);
			}
		}
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
		this.densities.fill(selfDensity);
		this.gradP.forEach(v => v.set(0, 0, 0));
		this.lapV.forEach(v => v.set(0, 0, 0));
		this.calcProperties();
		for (let i = 0; i < this.n; i++) {
			this.v[i].add(this.gravity(i).multiplyScalar(dt));
			this.v[i].add(this.gradP[i].clone().multiplyScalar(dt));
			this.v[i].add(this.lapV[i].clone().multiplyScalar(dt * this.viscosity));
			const prevPos = this.r[i].clone();
			this.r[i].add(this.v[i].clone().multiplyScalar(dt));
			if (this.r[i].y < worldBottom) {
				this.r[i].y = worldBottom;
				this.v[i].y = -this.v[i].y * restitution;
			}
			if (Math.abs(this.r[i].x) > worldWidth) {
				this.r[i].x = Math.sign(this.r[i].x) * (1 + restitution) * worldWidth - this.r[i].x * restitution;
				this.v[i].x *= -restitution;
				/* this.v[i].x = this.v[i].x * (0.04 ** dt) -
					Math.sign(this.r[i].x) * (Math.abs(this.r[i].x) - worldWidth) * dt; */
			}
			if (Math.abs(this.r[i].y) > worldWidth) {
				this.r[i].y = Math.sign(this.r[i].y) * (1 + restitution) * worldWidth - this.r[i].y * restitution;
				this.v[i].y *= -restitution;
			}
			if (Math.abs(this.r[i].z) > worldWidth) {
				this.r[i].z = Math.sign(this.r[i].z) * (1 + restitution) * worldWidth - this.r[i].z * restitution;
				this.v[i].z *= -restitution;
			}
			this.updateChunk(i, prevPos);
		}
	}
}

export const Simulation = new FluidHandler(2000);

// @ts-ignore
window.Simulation = Simulation;