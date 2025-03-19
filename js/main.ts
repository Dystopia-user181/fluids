import * as THREE from "three";
import { Simulation, worldBottom } from "./fluid.js";
import { scene } from "./scene.js";

import "./lights.js";

const camera = new THREE.PerspectiveCamera(90, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.z = 13;

const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);


const dropletRadius = 0.08;
const spheres: THREE.Mesh<THREE.SphereGeometry, THREE.MeshBasicMaterial>[] = [];
for (let i = 0; i < Simulation.n; i++) {
	const fluidShape = new THREE.SphereGeometry(dropletRadius, 16, 12);
	const fluidColour = new THREE.MeshBasicMaterial({ color: 0x0f659e });
	spheres.push(new THREE.Mesh(fluidShape, fluidColour));
	scene.add(spheres[spheres.length - 1]);
}

const geometry = new THREE.PlaneGeometry(10000, 10000);
const material = new THREE.MeshPhongMaterial({ color: 0x777777, side: THREE.DoubleSide });
const plane = new THREE.Mesh(geometry, material);
plane.position.set(0, worldBottom, 0);
plane.rotation.x += Math.PI / 2;
scene.add(plane);

const ticks = Array(30).fill(0.0166666) as number[];
let time = 0;
function animate(newTime: number) {
	const dt = (newTime - time) / 1000;
	ticks.shift();
	ticks.push(dt);
	time = newTime;
	Simulation.tick(dt);
	for (let i = 0; i < Simulation.n; i++) {
		const d = Simulation.densities[i] * 2 - 1;
		const hue = 204 / (1 + d * d * d * d * d * d);
		spheres[i].position.copy(Simulation.r[i]);
		spheres[i].material.color.setHSL(hue / 360, 0.83, 0.34);
	}
	renderer.render(scene, camera);
	(document.getElementById("fps-text") as HTMLDivElement).innerText =
		`fps: ${(30 / ticks.reduce((a, b) => a + b)).toFixed(2)}`;
}

window.addEventListener("keypress", () => renderer.setAnimationLoop(animate));