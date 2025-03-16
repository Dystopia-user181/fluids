import { resolve } from "path";

import { defineConfig } from "vite";

import checker from "vite-plugin-checker";
import eslint from "vite-plugin-eslint";
import tsconfigPaths from "vite-tsconfig-paths";

// https://vitejs.dev/config/
export default defineConfig({
	plugins: [
		eslint({
			failOnError: false
		}),
		checker({ typescript: true }),
		tsconfigPaths()
	],
	base: "./",
	resolve: {
		alias: {
			"@": resolve(__dirname, "./src")
		},
	},
});