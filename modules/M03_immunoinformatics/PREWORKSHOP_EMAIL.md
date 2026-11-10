# Pre workshop email, module M03.1 Immunoinformatics

Draft to send to participants. Adjust the dates and links before sending.

---

**Subject: Tutorial IP3, what to install before the immunoinformatics session**

Hello,

Module M03.1 of the tutorial goes from a protein sequence to a ranked shortlist
of T cell epitope candidates. We will work on **fhaC**, the filamentous
hemagglutinin transporter of *Bordetella pertussis*. The sequence is attached
and also ships with the repository, so you do not need to download it anywhere.

**You can follow the whole session with nothing installed.** There is a Colab
notebook that runs every core exercise in the browser. If that is your plan,
you only need a Google account and you can stop reading here.

If you want to run the real tools on your own machine, read on. Several of them
need an academic licence, and the requests can take a couple of days to come
back, so start now rather than the night before.

---

## If you want the local setup

**Base environment.** Python 3.10 or newer, then from the module folder:

```bash
pip install -r requirements.txt
```

```bash
mhcflurry-downloads fetch models_class1_presentation
```

```bash
python utils/download_data.py
```

Check it worked:

```bash
python solutions/ex11_topology_reference.py
```

It should print a signal peptide at residues 1 to 30 and 16 transmembrane beta
strands. If you see that, you are ready for most of the session.

**Licensed tools, request these early.**

- **NetMHCpan 4.1** and **NetMHCIIpan 4.0**, from DTU Health Tech. Academic
  licence, manual download.
- **SignalP 6.0**, also DTU, same process.
- **DeepTMHMM**, from DTU as well. Note that the BioLib route now needs an
  account and no longer runs locally.
- **MixMHCpred 2.2** and **PRIME 2.0**, from the Gfeller lab on GitHub. Clone
  the `v2.2` and `v2.0` tags, not the default branches.
- **DeepImmuno**, a plain GitHub clone, no licence needed.

Full instructions, including the environment variables each wrapper expects,
are in `QUICKSTART.md`.

**Linux is the smooth path** for the licensed tools. On macOS most of them work
with some effort. On Windows, use WSL or stay with Colab.

---

## What to bring

Nothing beyond a laptop. If the installs defeat you, come anyway. The Colab
notebook covers slides 18, 19, 22, 23 and 24 and the repository ships a
verified run of DeepTMHMM on fhaC, so no exercise is blocked by a missing tool.

Bring one thing in your head though. Before we run anything, we will ask you to
predict what each tool should say about fhaC. Being wrong on paper first is the
point of the session.

See you there,
The BigMind team
