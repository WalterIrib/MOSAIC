#!/bin/bash
# B1 — Sync Chimeric.out.junction + DCC output nel canonico
# Per 21 campioni (9 De Cecco + 12 METACAR).
#
# Cosa copia:
#   mapping/bam/<SID>/<SID>_Chimeric.out.junction   (junction file per pipeline M06)
#   circrna/<SID>/CircRNACount                        (DCC output)
#   circrna/<SID>/CircCoordinates                     (DCC output)
#   circrna/<SID>/dcc.log                             (DCC run log, per tracciabilità)
#
# Cosa NON copia:
#   *.circRNA, *.circRNAmapped  (file intermedi DCC, non fanno parte dell'output)
#   Chimeric.out.sam            (già presente nel canonico — MAI sovrascritto)
#
# MD5 check: Chimeric.out.sam pre/post deve essere identico.
# Lancia da: /home/walter/Desktop/mosaic1.1/
# Log: scripts/b1_sync.log

set -euo pipefail
LOG=/home/walter/Desktop/mosaic1.1/scripts/b1_sync.log
exec > >(tee -a "$LOG") 2>&1

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
ok()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✓ $*"; }
err() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ✗ $*"; }

log "=== B1 sync DCC → canonico ==="

# ─────────────────────────────────────────────────────────────────────────────
# Dataset definitions
# ─────────────────────────────────────────────────────────────────────────────

declare -A CANONICAL SRC_STAR SRC_DCC

# De Cecco (9 campioni)
DECECCA_CANONICAL=/home/walter/Desktop/DECECCA2019/mosaic_output
DECECCA_SRC_BASE=/home/walter/Desktop/DECECCA2019/mosaic_output/dcc
for SID in EP1 EP2 EP3 SEN_E1 SEN_E2 SEN_E3 SEN_L1 SEN_L2 SEN_L3; do
    CANONICAL[$SID]=$DECECCA_CANONICAL
    SRC_STAR[$SID]=$DECECCA_SRC_BASE/$SID/star
    SRC_DCC[$SID]=$DECECCA_SRC_BASE/$SID/dcc
done

# METACAR (12 campioni)
METACAR_CANONICAL=/home/walter/Desktop/METACAR/mosaic_output
METACAR_SRC_BASE=/home/walter/Desktop/METACAR/mosaic_output/dcc_rerun
for SID in ChaRNA1 ChaRNA2 ChaRNA3 \
            mito-ChaRNA1 mito-ChaRNA2 mito-ChaRNA3 \
            mito-totRNA1 mito-totRNA2 mito-totRNA3 \
            totRNA1 totRNA2 totRNA3; do
    CANONICAL[$SID]=$METACAR_CANONICAL
    SRC_STAR[$SID]=$METACAR_SRC_BASE/$SID/star
    SRC_DCC[$SID]=$METACAR_SRC_BASE/$SID/dcc
done

ALL_SAMPLES=(
    EP1 EP2 EP3 SEN_E1 SEN_E2 SEN_E3 SEN_L1 SEN_L2 SEN_L3
    ChaRNA1 ChaRNA2 ChaRNA3
    mito-ChaRNA1 mito-ChaRNA2 mito-ChaRNA3
    mito-totRNA1 mito-totRNA2 mito-totRNA3
    totRNA1 totRNA2 totRNA3
)

ERRORS=0

# ─────────────────────────────────────────────────────────────────────────────
# FASE 1: Pre-copy MD5 di tutti i Chimeric.out.sam
# ─────────────────────────────────────────────────────────────────────────────
log "--- FASE 1: MD5 pre-copy di Chimeric.out.sam ---"
declare -A SAM_MD5_PRE

for SID in "${ALL_SAMPLES[@]}"; do
    SAM=${CANONICAL[$SID]}/mapping/bam/$SID/${SID}_Chimeric.out.sam
    if [[ -f "$SAM" ]]; then
        SAM_MD5_PRE[$SID]=$(md5sum "$SAM" | awk '{print $1}')
        log "$SID SAM MD5 pre: ${SAM_MD5_PRE[$SID]}"
    else
        log "$SID: SAM non trovato in canonico — skip MD5 pre"
    fi
done

# ─────────────────────────────────────────────────────────────────────────────
# FASE 2: Controlla file intermedi nei sorgenti (mai copiare .circRNA/.circRNAmapped)
# ─────────────────────────────────────────────────────────────────────────────
log "--- FASE 2: verifica file intermedi DCC nei sorgenti ---"
for SID in "${ALL_SAMPLES[@]}"; do
    INTERMEDIATES=$(find "${SRC_STAR[$SID]}" \
        \( -name "*.circRNA" -o -name "*.circRNAmapped" \) 2>/dev/null)
    if [[ -n "$INTERMEDIATES" ]]; then
        log "$SID: file intermedi trovati (NON verranno copiati): $(echo $INTERMEDIATES | tr '\n' ' ')"
    fi
done
log "Verifica intermedi completata — nessuno sarà sincronizzato."

# ─────────────────────────────────────────────────────────────────────────────
# FASE 3: Sync per campione
# ─────────────────────────────────────────────────────────────────────────────
log "--- FASE 3: sync ---"
SYNCED=0

for SID in "${ALL_SAMPLES[@]}"; do
    log "──── $SID ────"
    CAN=${CANONICAL[$SID]}
    BAM_DIR=$CAN/mapping/bam/$SID
    CIRC_DIR=$CAN/circrna/$SID

    # Destinazioni
    JXN_DST=$BAM_DIR/${SID}_Chimeric.out.junction
    CNT_DST=$CIRC_DIR/CircRNACount
    CRD_DST=$CIRC_DIR/CircCoordinates
    LOG_DST=$CIRC_DIR/dcc.log

    # Sorgenti
    JXN_SRC=${SRC_STAR[$SID]}/${SID}_Chimeric.out.junction
    CNT_SRC=${SRC_DCC[$SID]}/CircRNACount
    CRD_SRC=${SRC_DCC[$SID]}/CircCoordinates
    LOG_SRC=${SRC_DCC[$SID]}/dcc.log

    # Verifica sorgenti
    MISSING=0
    for F in "$JXN_SRC" "$CNT_SRC" "$CRD_SRC"; do
        if [[ ! -f "$F" ]]; then
            err "$SID: sorgente mancante: $F"
            MISSING=1
        fi
    done
    if [[ $MISSING -eq 1 ]]; then
        err "$SID: SALTATO per sorgenti mancanti"
        ((ERRORS++)) || true
        continue
    fi

    # Verifica che junction non sia già presente nel canonico
    if [[ -f "$JXN_DST" && ! -L "$JXN_DST" ]]; then
        log "$SID: junction già presente nel canonico — verifico MD5"
        MD5_SRC=$(md5sum "$JXN_SRC" | awk '{print $1}')
        MD5_DST=$(md5sum "$JXN_DST" | awk '{print $1}')
        if [[ "$MD5_SRC" == "$MD5_DST" ]]; then
            ok "$SID: junction identico — skip copia"
        else
            err "$SID: junction nel canonico DIVERSO dalla sorgente — verifica manuale richiesta"
            ((ERRORS++)) || true
            continue
        fi
    fi

    # Crea destinazioni (mkdir -p crea bam/SID se non esiste — METACAR non l'ha)
    mkdir -p "$BAM_DIR" "$CIRC_DIR"

    # Copia junction (solo il file .junction — mai .circRNA/.circRNAmapped)
    if [[ ! -f "$JXN_DST" ]]; then
        cp "$JXN_SRC" "$JXN_DST"
        ok "$SID: junction copiato ($(du -sh "$JXN_DST" | cut -f1))"
    fi

    # Copia DCC output (CircRNACount, CircCoordinates, dcc.log)
    cp "$CNT_SRC" "$CNT_DST"
    cp "$CRD_SRC" "$CRD_DST"
    [[ -f "$LOG_SRC" ]] && cp "$LOG_SRC" "$LOG_DST" || true

    N_CIRC=$(awk 'NR>1 && $4>=2' "$CNT_DST" | wc -l)
    ok "$SID: CircRNACount ($N_CIRC circRNA >=2 jxn), CircCoordinates"

    ((SYNCED++)) || true
done

log ""
log "--- FASE 4: MD5 post-copy di Chimeric.out.sam ---"
log "(Solo campioni con SAM canonico pre-esistente — De Cecco. METACAR: SAM non in canonico, skip.)"
MISMATCH=0
for SID in "${ALL_SAMPLES[@]}"; do
    SAM=${CANONICAL[$SID]}/mapping/bam/$SID/${SID}_Chimeric.out.sam
    if [[ -n "${SAM_MD5_PRE[$SID]:-}" && -f "$SAM" ]]; then
        MD5_POST=$(md5sum "$SAM" | awk '{print $1}')
        if [[ "$MD5_POST" == "${SAM_MD5_PRE[$SID]}" ]]; then
            ok "$SID: Chimeric.out.sam invariato (MD5 OK)"
        else
            err "$SID: Chimeric.out.sam MD5 CAMBIATO — ANOMALIA CRITICA"
            ((MISMATCH++)) || true
        fi
    fi
done

log ""
log "=== B1 COMPLETATA ==="
log "Campioni sincronizzati: $SYNCED / ${#ALL_SAMPLES[@]}"
[[ $MISMATCH -gt 0 ]] && err "MD5 mismatch su $MISMATCH campioni — verifica manuale" || ok "MD5 check SAM: tutti invariati"
[[ $ERRORS -gt 0 ]]   && err "Errori: $ERRORS campioni con problemi" || ok "Zero errori"
