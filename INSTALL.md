# Installing LARMOR

A step-by-step guide to get LARMOR running on **Windows, macOS, or Linux**.
It takes about 10 minutes. If something goes wrong, jump to
[Troubleshooting](#troubleshooting) — it lists the errors people actually hit.

> **The one rule that avoids 90 % of problems:** use **Python 3.11** (3.10–3.12
> are fine; **not** 3.13 yet). One of LARMOR's dependencies, `mrsimulator`, is a
> compiled package that does not yet publish ready-made installers ("wheels") for
> the newest Python, so on 3.13 it tries to build itself from source and usually
> fails. The Conda route below pins the right Python for you automatically.

---

## Easiest — Windows one-click

If you're on Windows, you don't need to type any commands:

1. Get the code — `git clone https://github.com/sams808/LARMOR.git`, or download
   the ZIP from GitHub and unzip it.
2. Open the `LARMOR` folder and **double-click `install.bat`**.
3. A console window walks through the install (a few minutes). When it says
   *Done*, there's a **`LARMOR.bat` on your Desktop** — **double-click it to
   launch the app**.

`install.bat` uses your Conda if you have it (recommended), otherwise a Python
3.11 virtual environment. To update later, **double-click `update.bat`** (pulls
the latest code and refreshes everything). Both are safe to re-run.

> If `install.bat` closes instantly or Windows SmartScreen warns you: it's an
> unsigned local script. Click *More info → Run anyway*, or right-click
> `install.ps1 → Run with PowerShell`. If Conda isn't found and you have no
> Python 3.11, it will tell you to install Miniconda first (below).

The manual routes below do exactly the same thing by hand, and are what
macOS/Linux users follow.

---

## What you need first

- **~2 GB free disk** and an internet connection.
- Either **Miniconda** (recommended — handles the tricky compiled packages for
  you) **or** a **Python 3.11** installation with `pip`.
- **Git** (optional — you can also download the repository as a ZIP).

Get the code first, from a terminal:

```
git clone https://github.com/sams808/LARMOR.git
cd LARMOR
```

(Or download the ZIP from GitHub, unzip it, and open a terminal **inside** the
unzipped `LARMOR` folder.)

Everything below is run **from inside that `LARMOR` folder.**

---

## Option A — Conda (recommended, most reliable)

This is the surest route: Conda installs the compiled scientific packages as
prebuilt binaries, so nothing has to compile on your machine.

1. **Install Miniconda** if you don't have it:
   <https://docs.conda.io/en/latest/miniconda.html>
   (On Windows, this gives you an **"Anaconda Prompt"** in the Start menu — use
   that window for the next steps, *not* the regular Command Prompt.)

2. **Create the environment** (one time, from the `LARMOR` folder):

   ```
   conda env create -f environment.yml
   ```

   This makes an environment called **`larmor`** with the correct Python and
   every dependency, and installs LARMOR itself. It takes a few minutes.

3. **Activate it** whenever you want to use LARMOR:

   ```
   conda activate larmor
   ```

4. **Launch** (see [Launching](#launching) below).

If you already made the environment earlier and just pulled new code, refresh it
with `conda env update -f environment.yml`.

---

## Option B — Python + pip (no Conda)

Use this only if you already have **Python 3.11** (check with `python --version`).

1. **Create and activate a virtual environment** (from the `LARMOR` folder):

   - **Windows:**
     ```
     python -m venv .venv
     .venv\Scripts\activate
     ```
   - **macOS / Linux:**
     ```
     python3 -m venv .venv
     source .venv/bin/activate
     ```

2. **Install LARMOR and everything it needs:**

   ```
   pip install -r requirements.txt
   ```

   (Equivalently: `pip install -e ".[desktop]"`. The `[desktop]` part is what
   pulls in the graphical interface — don't leave it off.)

---

## Launching

Make sure the environment is **active** first (`conda activate larmor`, or your
`.venv` activated), then:

- **Any platform — from the terminal:**
  ```
  larmor desktop
  ```
- **Windows shortcut:** double-click **`LARMOR.bat`** in the `LARMOR` folder. It
  finds the `larmor` Conda environment automatically and starts the app.

The desktop window should open within a few seconds (the first launch is a little
slower while it warms up). You can also use LARMOR without the GUI:
`larmor info <path-to-data>`, `larmor fit <recipe.json>`, etc.

---

## Verify it worked

With the environment active:

```
larmor info --help
```

should print usage text (not an error). If the desktop window opens with
`larmor desktop`, you're done.

---

## Troubleshooting

**`conda: command not found` / `'conda' is not recognized`**
Conda isn't on your PATH. On Windows, open the **Anaconda Prompt** (Start menu)
instead of the plain Command Prompt. On macOS/Linux, close and reopen the
terminal after installing Miniconda, or run `source ~/miniconda3/bin/activate`.

**`No module named 'larmor'`**
The package isn't installed in the active environment. Make sure the environment
is active, then re-run the install command from inside the `LARMOR` folder
(`conda env create -f environment.yml`, or `pip install -r requirements.txt`).

**`No module named 'PySide6'` (or `pyqtgraph`) when running `larmor desktop`**
You installed only the core, not the graphical interface. Fix it with:
```
pip install -e ".[desktop]"
```
or, for a Conda env, `conda env update -f environment.yml`.

**`mrsimulator` fails to install / build**, e.g.
`error: Microsoft Visual C++ 14.0 or greater is required`, or
`Failed building wheel for mrsimulator`, or a long C/Cython compiler error:
this almost always means **your Python is too new** (3.13+) so pip couldn't find
a prebuilt wheel and tried to compile from source. Two fixes:
1. **Use Conda (Option A)** — it installs `mrsimulator` as a prebuilt binary and
   sidesteps compiling entirely. This is the recommended fix.
2. Or use **Python 3.11** for your virtual environment instead of the newest one.

**`LARMOR.bat` flashes open and closes**, or says
*"Could not find the 'larmor' conda environment"*: you haven't created the
environment yet, or it isn't named `larmor`. Do Option A step 2 first. (Run the
`.bat` from a terminal to read the message: `.\LARMOR.bat`.)

**Linux: the window doesn't appear / `xcb` or `libGL` errors.** Install the Qt
runtime libraries, e.g. on Debian/Ubuntu:
```
sudo apt install libgl1 libegl1 libxcb-cursor0
```

**The app opens but text is invisible / white-on-white.** Update to the latest
code (`git pull` then `conda env update -f environment.yml`) — LARMOR forces a
readable light theme, and this was fixed in an earlier version.

**Still stuck?** Run `larmor desktop` from the terminal and copy the full error
message — it names the missing piece. Also confirm your Python version with
`python --version` (aim for 3.11).

---

*LARMOR — Sam Soudani, McCloy group, Washington State University.*
