import fs from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const repoRoot = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const sourceRoot = path.join(repoRoot, "runtime");
const targetRoot = path.join(repoRoot, "src-tauri", "resources", "runtime");
const ignoredDirectories = new Set([".nextflow", "__pycache__", ".pytest_cache"]);

function shouldCopy(name, isDirectory) {
  if (isDirectory && ignoredDirectories.has(name)) return false;
  if (/^\.nextflow\.log(?:\.\d+)?$/.test(name)) return false;
  if (/\.(?:pyc|pyo)$/.test(name)) return false;
  if (name === ".DS_Store") return false;
  return true;
}

async function copyDirectory(source, target) {
  await fs.mkdir(target, { recursive: true });
  for (const entry of await fs.readdir(source, { withFileTypes: true })) {
    if (!shouldCopy(entry.name, entry.isDirectory())) continue;
    const sourcePath = path.join(source, entry.name);
    const targetPath = path.join(target, entry.name);
    if (entry.isDirectory()) {
      await copyDirectory(sourcePath, targetPath);
    } else if (entry.isFile()) {
      await fs.copyFile(sourcePath, targetPath);
    }
  }
}

await fs.rm(targetRoot, { recursive: true, force: true });
await copyDirectory(sourceRoot, targetRoot);
await fs.writeFile(path.join(targetRoot, ".gitkeep"), "\n", "utf8");

const copied = await fs.readdir(targetRoot);
if (!copied.includes("main.nf") || !copied.includes("environment.yml")) {
  throw new Error("Desktop runtime preparation did not copy the required workflow files.");
}

console.log(`Prepared desktop runtime resource at ${targetRoot}`);
