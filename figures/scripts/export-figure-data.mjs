#!/usr/bin/env node

import { mkdirSync, readFileSync, writeFileSync } from "node:fs";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";
import { spawnSync } from "node:child_process";

const rootDir = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const envPath = resolve(rootDir, ".env");
const sqlPath = resolve(rootDir, "sql", "figure_data.sql");
const outputPath = resolve(rootDir, "data", "figures.json");

function readDotEnv(path) {
  const values = {};
  const text = readFileSync(path, "utf8");

  for (const line of text.split(/\r?\n/)) {
    const trimmed = line.trim();
    if (!trimmed || trimmed.startsWith("#")) continue;

    const match = trimmed.match(/^([A-Za-z_][A-Za-z0-9_]*)=(.*)$/);
    if (!match) continue;

    let value = match[2].trim();
    if (
      (value.startsWith("\"") && value.endsWith("\"")) ||
      (value.startsWith("'") && value.endsWith("'"))
    ) {
      value = value.slice(1, -1);
    }

    values[match[1]] = value;
  }

  return values;
}

function redact(text, secret) {
  if (!secret) return text;
  return text.split(secret).join("[DATABASE_URL]");
}

const dotEnv = readDotEnv(envPath);
const databaseUrl = process.env.DATABASE_URL || dotEnv.DATABASE_URL;

if (!databaseUrl) {
  console.error("DATABASE_URL was not found in .env or the process environment.");
  process.exit(1);
}

const result = spawnSync(
  "psql",
  [
    databaseUrl,
    "-v",
    "ON_ERROR_STOP=1",
    "-X",
    "-q",
    "-t",
    "-A",
    "-f",
    sqlPath
  ],
  {
    cwd: rootDir,
    encoding: "utf8",
    env: {
      ...process.env,
      ...dotEnv,
      PGCONNECT_TIMEOUT: process.env.PGCONNECT_TIMEOUT || "15"
    },
    maxBuffer: 64 * 1024 * 1024
  }
);

if (result.status !== 0) {
  const stderr = redact(result.stderr || "", databaseUrl);
  const stdout = redact(result.stdout || "", databaseUrl);
  console.error(stderr || stdout || "psql failed without output.");
  process.exit(result.status || 1);
}

const jsonText = result.stdout.trim();
let parsed;

try {
  parsed = JSON.parse(jsonText);
} catch (error) {
  console.error("The SQL export did not produce valid JSON.");
  console.error(error.message);
  process.exit(1);
}

mkdirSync(dirname(outputPath), { recursive: true });
writeFileSync(outputPath, `${JSON.stringify(parsed, null, 2)}\n`);

console.log(`Wrote ${outputPath}`);
console.log(
  [
    `${parsed.totals.entity_count.toLocaleString()} entities`,
    `${parsed.totals.relation_count.toLocaleString()} relations`,
    `${parsed.totals.built_resource_count.toLocaleString()} built resources`
  ].join(" | ")
);
