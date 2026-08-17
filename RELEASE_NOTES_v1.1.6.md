## v1.1.6 — Fix: correct cfg attribute in provenance logging

**Tag:** v1.1.6 — commit `372ef2f`
**Date:** 2026-08-17
**Supersedes:** v1.1.5 (ritirato)

### Bug fixed

v1.1.5 (cd29e9f) aggiungeva una riga `[PROVENANCE]` come prima riga di
`mosaic_output/logs/session.log` a ogni `mosaic run`, ma usava
`cfg.wd` (attributo inesistente) invece di `cfg.project.wd`.

**Errore:** `AttributeError: 'MosaicConfig' object has no attribute 'wd'`

v1.1.6 corregge l'attributo. La funzionalità di provenance logging è
ora operativa: la prima riga di `session.log` è:

```
[PROVENANCE] <ISO-timestamp>
  mosaic.__file__=<path>/mosaic_code/mosaic/cli.py
  scripts_dir=<path>/mosaic_code/scripts
```

### Perché questa funzionalità

La provenienza del codice dipendeva dalla directory di lancio (`cwd` in
`sys.path[0]`): se il runner veniva eseguito da una directory contenente
una cartella `mosaic/`, Python trovava quella cartella prima del
pacchetto installato. Loggare `mosaic.__file__` e `scripts_dir` a
runtime rende la provenienza verificabile a posteriori senza ricostruire
lo stato dell'ambiente.

### Files changed

- `mosaic/cli.py`: `cfg.wd` → `cfg.project.wd` in `cmd_run`

### Does this change pipeline results?

No. La riga di provenance viene scritta prima che qualsiasi modulo
venga avviato. Nessun risultato biologico è influenzato.

### Tag e citazione

Il manoscritto cita **v1.1.2** (commit `e022cf6`), invariato.
v1.1.6 non è citato nel manoscritto.

### Nota su smoke test

Il crash di v1.1.5 ha evidenziato l'assenza di un test end-to-end che
verifichi l'avvio della pipeline su un config minimale. Un test di
questo tipo (`mosaic run` su 2 campioni finti, 1 confronto, moduli
ridotti al minimo) sarà aggiunto in una release successiva come
`tests/test_smoke_run.py`.
