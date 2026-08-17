## v1.1.5 — RETIRED (non avviabile)

**Tag:** v1.1.5 — commit `cd29e9f`
**Status:** RITIRATO — tag eliminato da GitHub il 2026-08-17
**Sostituito da:** v1.1.6 (372ef2f)

### Perché ritirato

v1.1.5 introduceva il logging di provenance in `cmd_run` (cli.py), ma
usava `cfg.wd` invece di `cfg.project.wd`. L'attributo `wd` non esiste
in `MosaicConfig` — si trova in `MosaicConfig.project.wd`.

Effetto: `AttributeError: 'MosaicConfig' object has no attribute 'wd'`
a ogni avvio di `mosaic run`, su qualunque config. La pipeline non
arrivava al primo modulo.

v1.1.5 è stato pushato e immediatamente ritirato lo stesso giorno
(2026-08-17) senza mai essere usato in alcun run. Il commit `cd29e9f`
rimane nella storia (è il padre di v1.1.6) ma il tag è eliminato.

### Fix

v1.1.6 (372ef2f) corregge l'attributo in `cfg.project.wd`.
