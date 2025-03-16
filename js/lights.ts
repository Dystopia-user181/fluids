import * as THREE from "three";
import { scene } from "./scene.js";

const intensity = 500;
const light = new THREE.PointLight(0xffffff, intensity);
scene.add(light);
light.position.set(0, 10, 0);
const ambientLight = new THREE.AmbientLight(0x707070);
scene.add(ambientLight);