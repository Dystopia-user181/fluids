import { Vector3 } from "three";

export const worldBottom = -16;
const worldWidth = 3;
const influence = 0.5;
const influenceSt = 1.5;
const influenceSq = influence * influence;
const influenceStSq = influenceSt * influenceSt;
const dimension = 2;

const pressureForceStrength = 5e6;
const constantEncodeOffset = 10;

function interLeave3(input: number) {
	let x = input;
	x = (x | x << 8) & 0xf00f00f00f00f;
	x = (x | x << 4) & 0x30c30c30c30c3;
	x = (x | x << 2) & 0x0249249249249;
	return x;
}
function getChunkWKey(pos: Vector3) {
	return (interLeave3(Math.round((pos.x + constantEncodeOffset) / influence)) +
		(interLeave3(Math.round((pos.y + constantEncodeOffset) / influence)) << 1) +
		(interLeave3(Math.round((pos.z + constantEncodeOffset) / influence)) << 2)) & 4095;
}
function getChunkWStKey(pos: Vector3) {
	return (interLeave3(Math.round((pos.x + constantEncodeOffset) / influenceSt)) +
		(interLeave3(Math.round((pos.y + constantEncodeOffset) / influenceSt)) << 1) +
		(interLeave3(Math.round((pos.z + constantEncodeOffset) / influenceSt)) << 2)) & 4095;
}

function getPressure(_density: number) {
	const density = _density / targetDensity - 0.85;
	const u1 = density * density * density;
	return u1 * u1 * density * pressureForceStrength;
}

const W = (() => {
	const normalize = 1 / (influence ** dimension);
	const invInfSq = 1 / influenceSq;
	return (magDrSq: number) => normalize * (magDrSq * invInfSq - 2 * Math.sqrt(magDrSq * invInfSq) + 1);
})();

const gradW = (() => {
	const normalize = 2 / (influence ** dimension);
	const invInfSq = 1 / influenceSq;
	return (dr: Vector3) => dr.clone().multiplyScalar(normalize * (invInfSq - Math.sqrt(invInfSq / dr.lengthSq())));
})();

const WSt = (() => {
	const normalize = 1 / (influenceSt ** dimension);
	const invInfSq = 1 / influenceStSq;
	const halfInf = influenceSt / 2;
	const halfInfSq = halfInf * halfInf;
	return (magDrSq: number) => normalize * (magDrSq < halfInfSq ? Math.sqrt(magDrSq * invInfSq)
		: magDrSq * invInfSq - 2 * Math.sqrt(magDrSq * invInfSq) + 1);
})();

const targetDensity = 1.5 * W(0), selfDensity = W(0);

interface BoundaryConfig {
	displacement: (x: Vector3) => Vector3,
	sgnDistance: (x: Vector3) => number,
	restitution: number,
	adhesion: number,
}
class Boundary {
	config: BoundaryConfig;

	constructor(x: Partial<BoundaryConfig>) {
		if (!x.displacement) throw "Boundary needs displacement function";
		if (!x.sgnDistance) throw "Boundary needs distance function";
		this.config = {
			displacement: x.displacement,
			sgnDistance: x.sgnDistance,
			restitution: x.restitution || 0,
			adhesion: x.adhesion || 0,
		};
	}

	distance(x: Vector3) {
		return this.config.sgnDistance(x.clone());
	}

	displacement(x: Vector3) {
		return this.config.displacement(x.clone());
	}

	handleCollision(r: Vector3, v: Vector3) {
		if (this.distance(r) >= 0) return;
		const d = this.displacement(r);
		r.sub(d);
		v.sub(v.clone().projectOnVector(d).multiplyScalar(1 + this.config.restitution));
	}
}

class FluidHandler {
	n = 0;
	r: Vector3[] = [];
	v: Vector3[] = [];
	chunkedR = Array(4096).fill(0).map(() => new Set<number>());
	chunkedRSt = Array(4096).fill(0).map(() => new Set<number>());
	inChunk: number[] = [];
	inChunkSt: number[] = [];

	boundaries: Boundary[] = [];

	densities: number[] = [];
	gradP: Vector3[] = [];
	lapV: Vector3[] = [];

	viscosity = 1;
	tension = 10;

	timeElapsed = 0;

	constructor(number: number) {
		this.n = number;
		for (let i = 0; i < this.n; i++) {
			this.r.push(new Vector3(Math.random() * 6 - 3, Math.random() * 32 - 3));
			this.v.push(new Vector3(0, 0, 0));
			this.densities.push(selfDensity);
			this.gradP.push(new Vector3(0, 0, 0));
			this.lapV.push(new Vector3(0, 0, 0));
			this.inChunk.push(4098);
			this.inChunkSt.push(4098);
			this.updateChunk(i);
		}
	}

	_updateChunk(keyFunc: (x: Vector3) => number, chunkList: Set<number>[], chunkIdentifier: number[],
		radius: number, u: number, _prevPos?: Vector3) {
		const newKey = keyFunc(this.r[u]);
		if (chunkIdentifier[u] === newKey && !_prevPos) return;
		const prevPos = _prevPos ? _prevPos : new Vector3(Infinity, Infinity, Infinity);
		const dx = Math.round((this.r[u].x + constantEncodeOffset) / radius) -
			Math.round((prevPos.x + constantEncodeOffset) / radius);
		const dy = Math.round((this.r[u].y + constantEncodeOffset) / radius) -
			Math.round((prevPos.y + constantEncodeOffset) / radius);
		const dz = Math.round((this.r[u].z + constantEncodeOffset) / radius) -
			Math.round((prevPos.z + constantEncodeOffset) / radius);
		if (_prevPos) {
			for (let i = -1; i < 2; i++) {
				for (let j = -1; j < 2; j++) {
					for (let k = -1; k < 2; k++) {
						// eslint-disable-next-line max-depth
						if (-dx + i >= -1 && -dx + i <= 1 &&
							-dy + j >= -1 && -dy + j <= 1 &&
							-dz + k >= -1 && -dz + k <= 1) continue;
						const newPos = new Vector3(
							prevPos.x + i * radius,
							prevPos.y + j * radius,
							prevPos.z + k * radius
						);
						const key = keyFunc(newPos);
						chunkList[key].delete(u);
					}
				}
			}
		}
		chunkIdentifier[u] = newKey;
		for (let i = -1; i < 2; i++) {
			for (let j = -1; j < 2; j++) {
				for (let k = -1; k < 2; k++) {
					// eslint-disable-next-line max-depth
					if (dx + i >= -1 && dx + i <= 1 &&
						dy + j >= -1 && dy + j <= 1 &&
						dz + k >= -1 && dz + k <= 1) continue;
					const newPos = new Vector3(
						this.r[u].x + i * radius,
						this.r[u].y + j * radius,
						this.r[u].z + k * radius
					);
					const key = keyFunc(newPos);
					chunkList[key].add(u);
				}
			}
		}
	}

	updateChunk(u: number, _prevPos?: Vector3) {
		this._updateChunk(getChunkWKey, this.chunkedR, this.inChunk, influence, u, _prevPos);
		this._updateChunk(getChunkWStKey, this.chunkedRSt, this.inChunkSt, influenceSt, u, _prevPos);
	}

	gravity(i: number) {
		// const lengthSq = this.r[i].x * this.r[i].x + 0.01 * (this.r[i].y + worldWidth) * (this.r[i].y + worldWidth);
		const phase = this.timeElapsed * 0.15;
		const mulFactor = Math.sin(phase) > 0.93 ? 0 : 0;
		const convection = new Vector3(mulFactor, 0, 0);
		return new Vector3(0, -3, 0).add(convection);
		const rr = this.r[i].clone().add(new Vector3(0, 0, 0));
		return rr.clone().multiplyScalar(-Math.min(3 / (rr.lengthSq() ** 1.5), 1));
		// const phase = this.timeElapsed * 0.5;
		const r = this.r[i].clone().add(new Vector3(4 * Math.cos(phase), 4 * Math.sin(phase), 0));
		const g1 = r.clone().multiplyScalar(-Math.min(0.1 / (r.lengthSq() ** 1.5), 1));
		const r2 = r.sub(new Vector3(4 * Math.cos(phase), 4 * Math.sin(phase), 0));
		const g2 = r2.clone().multiplyScalar(-Math.min(10 / (r2.lengthSq() ** 1.5), 1));
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
					if (magDrSq > influenceSq || magDrSq < 1e-15) continue;
					const dr = new Vector3(
						drx,
						dry,
						drz,
					);
					const density = W(magDrSq);
					this.densities[u] += density;
					forcesMag.push(gradW(dr));
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
				const multiplier = Math.min(getPressure(this.densities[u]), 1000);
				forcesMag[i].multiplyScalar(multiplier / this.densities[u] / this.densities[u]);
				this.gradP[u].add(forcesMag[i]);
				this.gradP[forcesBy[i]].sub(forcesMag[i]);
			}
		}
	}

	calcTension() {
		for (let i = 0; i < this.n; i++) {
			for (const boundary of this.boundaries) {
				const d = boundary.distance(this.r[i]);
				if (d < influenceSt) {
					this.gradP[i].sub(boundary.displacement(this.r[i]).normalize().multiplyScalar(
						WSt(d) * boundary.config.adhesion
					));
				}
			}
		}
		const sortedIds = Array.from({ length: this.n }, (_, i) => i)
			.sort((a, b) => this.inChunkSt[a] - this.inChunkSt[b]);
		const corrChunk = sortedIds.map(x => this.inChunkSt[x]);
		let ptr = 0;
		while (ptr < this.n) {
			const chunkBegins = ptr;
			const chunk = corrChunk[ptr];
			while (corrChunk[ptr] === corrChunk[ptr + 1]) ptr++;
			ptr++;
			const chunkEnds = ptr;
			// This is to help with caching
			const corrR = sortedIds.slice(chunkBegins, chunkEnds).map(x => this.r[x].clone());
			const tensions = corrR.map(() => new Vector3());
			for (const v of this.chunkedRSt[chunk]) {
				const { x: r1x, y: r1y, z: r1z } = this.r[v];
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
					if (magDrSq > influenceStSq || magDrSq < 1e-15) continue;
					const dr = new Vector3(
						drx,
						dry,
						drz,
					);
					const densitySt = WSt(magDrSq);
					tensions[i - chunkBegins].add(dr.normalize().multiplyScalar(
						Math.min(densitySt * this.tension, 100)
					));
				}
			}
			for (let i = chunkBegins; i < chunkEnds; i++) {
				const u = sortedIds[i];
				this.gradP[u].add(tensions[i - chunkBegins]);
			}
		}
	}

	tick(dt: number) {
		const _dt = Math.min(dt, 0.03);
		const n = Math.ceil(_dt / 0.01);

		for (let i = 0; i < n; i++) {
			// let prev = calledTimes;
			this._tick(_dt / n);
			// console.log(calledTimes - prev)
		}
	}

	_tick(dt: number) {
		this.timeElapsed += dt;
		this.densities.fill(selfDensity);
		this.gradP.forEach(v => v.set(0, 0, 0));
		this.lapV.forEach(v => v.set(0, 0, 0));
		this.calcProperties();
		this.calcTension();
		for (let i = 0; i < this.n; i++) {
			this.v[i].add(this.gravity(i).multiplyScalar(dt));
			this.v[i].add(this.gradP[i].clone().multiplyScalar(dt));
			this.v[i].add(this.lapV[i].clone().multiplyScalar(dt * this.viscosity));
			const prevPos = this.r[i].clone();
			this.r[i].add(this.v[i].clone().multiplyScalar(dt));
			for (const boundary of this.boundaries) {
				boundary.handleCollision(this.r[i], this.v[i]);
			}
			this.updateChunk(i, prevPos);
		}
	}
}

export const Simulation = new FluidHandler(1000);
Simulation.boundaries.push(new Boundary({
	sgnDistance(r: Vector3) { return worldWidth - r.x; },
	displacement(r: Vector3) { return new Vector3(r.x - worldWidth, 0, 0); },
	restitution: 0.1,
	adhesion: 100,
}),
new Boundary({
	sgnDistance(r: Vector3) { return r.x + worldWidth; },
	displacement(r: Vector3) { return new Vector3(r.x + worldWidth, 0, 0); },
	restitution: 0.1,
	adhesion: 100,
}),
new Boundary({
	sgnDistance(r: Vector3) { return r.y + worldWidth; },
	displacement(r: Vector3) { return new Vector3(0, r.y + worldWidth, 0); },
	restitution: 0.1,
	adhesion: 1,
}));

// @ts-ignore
window.Simulation = Simulation;