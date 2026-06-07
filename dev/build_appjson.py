#!/usr/bin/env python
"""Regenerate app.json (the shinylive payload) from app.R.

The shinylive runtime in ./shinylive reads app.json, a JSON array of files:
[{"name": "app.R", "content": "<source>", "type": "text"}]

Run from the repo root:  python dev/build_appjson.py
"""
import json
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
app_r = os.path.join(ROOT, "app.R")
out = os.path.join(ROOT, "app.json")

with open(app_r, "r", encoding="utf-8") as f:
    content = f.read()

payload = [{"name": "app.R", "content": content, "type": "text"}]

with open(out, "w", encoding="utf-8", newline="\n") as f:
    json.dump(payload, f, ensure_ascii=False, indent=2)
    f.write("\n")

print(f"Wrote {out} ({len(content)} chars of app.R)")
