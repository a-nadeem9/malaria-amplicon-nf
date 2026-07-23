import { defineConfig } from "vite";
import { fileURLToPath } from "node:url";

const projectRoot = fileURLToPath(new URL(".", import.meta.url));

export default defineConfig({
  root: projectRoot,
  server: {
    host: "127.0.0.1",
    port: 1420,
    strictPort: true,
  },
});

