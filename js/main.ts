import * as THREE from "three";
import { Simulation, worldBottom } from "./fluid.js";
import { scene } from "./scene.js";

import "./lights.js";

const camera = new THREE.PerspectiveCamera(90, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.position.z = 5;

const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);


const dropletRadius = 0.08;
const fluidShape = new THREE.SphereGeometry(dropletRadius, 16, 12);
const fluidColour = new THREE.MeshBasicMaterial({ color: 0x0f659e });
const spheres = Simulation.r.map(() => new THREE.Mesh(fluidShape, fluidColour));
for (const s of spheres) {
	scene.add(s);
}

const geometry = new THREE.PlaneGeometry(10000, 10000);
const material = new THREE.MeshPhongMaterial({ color: 0x777777, side: THREE.DoubleSide });
const plane = new THREE.Mesh(geometry, material);
plane.position.set(0, worldBottom, 0);
plane.rotation.x += Math.PI / 2;
scene.add(plane);

let time = 0;
function animate(newTime: number) {
	const dt = (newTime - time) / 1000;
	time = newTime;
	Simulation.tick(dt);
	for (let i = 0; i < Simulation.n; i++) {
		spheres[i].position.copy(Simulation.r[i]);
	}
	renderer.render(scene, camera);
}
renderer.setAnimationLoop(animate);