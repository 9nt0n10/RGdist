# Report: limiti di SPATA2, FDR a zero, `rm_zero_infl` e arricchimenti come benchmark

Informazioni ricavate in questa sessione dall'analisi dei risultati del benchmark
(6 dataset Visium: heme_1000, stab_inj1, stab_inj2, gbm, perio, head_and_neck),
dal paper di SPATA2 (Kueckelhaus et al., *Nature Communications* 2024) e dal
tutorial ufficiale (themilolab.github.io/SPATA2). Dove un'informazione è una
deduzione e non un dato, è indicato.

**Pipeline**
- `pipeline_v8`: SPATA2 con `rm_zero_infl = TRUE`, output in `benchmark_v8/<dataset>`.
- `pipeline_v8_noinfl`: SPATA2 con `rm_zero_infl = FALSE`, output in
  `benchmark_v8_noinfl/<dataset>`.
- Per il resto le due pipeline sono identiche, e in entrambe ogni metodo entra con
  la propria lista di geni. Il filtro "espresso in ≥ 3 spot" vale solo per i GLM:
  SPATA2 usa il suo output intero, SpatialDE il suo CSV, Moran e Laplacian tutti i
  geni con almeno un conteggio.
- I numeri marcati "universo comune" vengono dalla versione precedente
  (`pipeline_v7`), dove tutti i metodi erano ristretti ai geni del filtro dei GLM e
  i geni non testati andavano in coda.

---

## 1. Come SPATA2 calcola il p-value (dal paper)

- **Statistica.** Per ogni gene, SPATA2 ricava il profilo di espressione in funzione
  della distanza dall'annotazione, lo smussa (loess) e ne misura la variazione
  totale (`tot_var`, TV). Più il profilo è liscio, più è basso il TV, e più per
  SPATA2 il gradiente è "non casuale".
- **Nullo.** "We first simulate random expression gradients by assigning random
  expression values to all data points within the screened area". La simulazione è
  ripetuta 10000 volte.
- **p-value.** È la frazione di gradienti simulati con TV minore o uguale a quello
  osservato. Prima si tolgono gli outlier del TV simulato con il criterio
  dell'intervallo interquartile (IQR). Poi si applica BH, e il risultato è la
  colonna `fdr`.
- **Conseguenza.** Il nullo è fatto di valori casuali **senza struttura spaziale**.
  Qualunque gene con una struttura spaziale, anche non legata all'annotazione, ha
  un profilo più liscio di tutti i profili casuali. Il p-value quindi non misura
  "gradiente legato a questa zona" ma "il profilo non è rumore puro". È lo stesso
  limite dei p-value parametrici del GLM.
- **Nel paper stesso,** un esempio sul glioblastoma dà 11.478 geni con FDR < 0.05.
- **Deduzione da verificare:** se il nullo è uguale per tutti i geni, il p-value è
  solo una trasformazione monotona di `tot_var` e non aggiunge informazione
  all'ordine.

## 2. Geni con p = 0 e FDR = 0

- **heme_1000, senza filtro:**
  - 1822 geni hanno p esattamente 0, quindi FDR 0, tutti pari merito;
  - 2422 geni hanno p ≤ 0.001.
- **gbm, `pipeline_v8_noinfl`:**
  - 3950 geni hanno p = 0;
  - le prime 3000 posizioni della classifica hanno tutte un solo valore di FDR (0);
  - fra 3000 e 4000 ci sono altri 950 geni con p = 0;
  - il primo FDR diverso da zero è 0.00048.
- **Tutto il blocco a FDR 0 è ordinato solo dallo spareggio,** cioè `tot_var`
  crescente (`arrange(fdr, tot_var)`).
- **Più simulazioni non risolvono.** Questi geni non stanno al limite della
  distribuzione nulla, stanno ben oltre. Con 100000 simulazioni avrebbero comunque
  p = 0.

## 3. Artefatto sui geni con molti zeri

- **heme_1000, senza filtro:**
  - l'output di SPATA2 ha 17941 geni;
  - 1747 di questi sono fuori dal filtro dei GLM (espressi in < 3 spot), e SPATA2 ne
    chiama significativo il **96.2%** (1681), contro l'11.3% dei geni espressi in ≥ 3
    spot;
  - fra i primi in classifica ci sono recettori olfattivi e vomeronasali
    (*Olfr718-ps1*, *Vmn1r32*), geni predetti (*Gm44928*) e geni testicolari (*Tex45*),
    senza senso biologico per un gradiente spaziale.
- **gbm, `pipeline_v8_noinfl`, composizione della classifica di SPATA2 per fasce:**

  | Rank | Geni espressi in < 3 spot | Geni del GS | `tot_var` mediana |
  |---|---|---|---|
  | 1–500 | 165 | 216 | 1.254 |
  | 501–1000 | 63 | 348 | 1.753 |
  | 1001–1500 | 319 | 87 | 2.000 |
  | 1501–2000 | 466 | 0 | 2.000 |
  | 2001–2500 | 425 | 15 | 2.000 |
  | 2501–3000 | 134 | 227 | 2.104 |

  La fascia 1500–2500 è il "buco" senza geni del Golden Standard nel boxplot dei
  rank: 891 geni su 1000 sono espressi in meno di 3 spot, tutti con p = 0 e
  `tot_var` = 2.000.
- **Spiegazione (dedotta).** SPATA2 riscala ogni profilo fra 0 e 1. Un gene acceso in
  un solo spot dà un profilo piatto con una punta che sale da 0 a 1 e riscende: la
  variazione totale è esattamente 2. I profili casuali sono più frastagliati, quindi
  il gene ha p = 0. Inoltre, nel blocco a FDR 0 lo spareggio per `tot_var` crescente
  premia proprio i profili piatti.
- **Gli autori stessi riconoscono il limite:** "Given the limitations of the
  spatial gradient screening algorithm in handling variables with a high proportion
  of zero values, we advise pre-emptively removing such variables."

## 4. `rm_zero_infl`: cosa fa

- **Regola (dal paper).** Per ogni gene si tolgono gli spot oltre 1.5 IQR sopra il
  terzo quartile o sotto il primo. Se restano solo zeri (tutti i valori non nulli
  sono outlier), il gene è "zero-inflated" e viene escluso dallo screening. Per quei
  geni non viene calcolato nessun p-value.
- **Effetto pratico (dedotto dalla regola).** Se un gene è zero in più del 75% degli
  spot, primo e terzo quartile valgono 0, l'IQR è 0 e ogni valore positivo è un
  outlier. Quindi l'opzione esclude in pratica i geni rilevati in meno di circa un
  quarto degli spot, cioè anche i geni espressi solo in una zona localizzata come una
  lesione o un aggregato cellulare.
- **Valore predefinito.** Il paper dice che l'opzione "can be activated by setting
  the parameter rm_zero_infl to TRUE", il che indica che il predefinito sia FALSE
  (deduzione). Il commento nella pipeline ("default del tool") è probabilmente
  sbagliato. Il tutorial ufficiale non nomina l'opzione.

## 5. SPATA2 con e senza `rm_zero_infl`

### gbm (`pipeline_v8` contro `pipeline_v8_noinfl`)

| | senza filtro | con filtro |
|---|---|---|
| Geni nella lista di SPATA2 | 20503 | 5455 (−73%) |
| Geni significativi (FDR < 0.05) | 8124 | 2419 |
| Geni del GS nella lista (su 7280) | 7280 | 3549 |
| GS recuperati (FDR < 0.05) | 3693 (50.7%) | 2199 (30.2%) |
| AUC sul GS, sulla propria lista | 0.633 | 0.847 |
| AUC sul GS, con i GS non testati in fondo | 0.633 | 0.413 |
| GLM lineare, AUC sul GS (per confronto) | 0.949 | 0.949 |

- **Geni del GS tolti dal filtro:** 3731. Fra questi ci sono marcatori di
  proliferazione (*Mki67*, *Cdk1*, *Plk1*, *Bub1*, *Cdc6*, *Ska1*, *Exo1*), di
  risposta all'interferone (*Apol9a*, *Apol9b*, *Sp140*, *Epsti1*) e altri geni del
  contesto tumorale (*Cdkn2a*, *Hmga2*, *Timp1*).
  - 1545 di questi erano significativi senza filtro, molti fra le prime 300
    posizioni.
  - Rank relativo mediano, senza filtro: 0.454 per i GS tolti, 0.332 per quelli tenuti.
- **L'ordine dei geni rimasti è identico** (Spearman 1): il filtro toglie geni, non
  cambia la classifica.
- **Sovrapposizione dei significativi fra i due run:** 2355 in comune, Jaccard 0.288.
- Fra i due run non cambia nessun altro metodo (controllo di coerenza superato).

### perio (`pipeline_v8`, con filtro)
- La lista di SPATA2 ha 8999 geni, contro i 15163 dei GLM, e solo 11 sono fuori dal
  filtro dei GLM. Nella versione con universo comune e senza filtro, SPATA2 aveva un
  FDR per 15149 geni.
- Toglie 126 dei 353 geni del Golden Standard, quasi tutti dell'infiltrato
  immunitario:
  - linfociti T: *CD3D*, *LCK*, *ZAP70*, *TRAC*, *TRBC1*, *CD247*, *IL7R*, *TCF7*,
    *CTLA4*;
  - linfociti B e plasmacellule: *CD79A*, *MZB1*, *IGHM*, *IGHG3*, *IRF4*, *POU2AF1*;
  - chemochine e marcatori di aggregati linfoidi: *CXCR4*, *CCR7*, *CCL19*, *CCL5*,
    *SELL*, *PTPRC*.

  Nella parodontite sono i geni dell'infiammazione, cioè il contesto biologico.

| AUC sul Golden Standard, perio | SPATA2 | GLM lineare | GLM con centro |
|---|---|---|---|
| sulla propria lista | 0.861 | 0.895 | 0.895–0.897 |
| solo sugli 8999 geni testati da SPATA2 | 0.861 | 0.915 | 0.915–0.927 |
| con i GS non testati in fondo | 0.554 | 0.895 | 0.895–0.897 |

### Copertura del Golden Standard senza filtro (gbm)
SPATA2 ha 7280 geni del GS nella lista, SpatialDE 7275, GLM, Moran e Laplacian
7267. Senza filtro SPATA2 calcola un p-value per praticamente tutti i geni rilevati
(20503 contro 20558 di Moran e Laplacian).

### In sintesi
- **Senza filtro:** SPATA2 testa tutti i geni, ma i geni con molti zeri ricevono
  p = 0 artefattuali e occupano posizioni alte (§3).
- **Con il filtro:** l'artefatto sparisce, ma insieme a esso sparisce una gran parte
  dei geni, compresi quelli localizzati del contesto biologico. L'AUC sulla
  propria lista sale solo perché quei geni del GS non vengono più valutati;
  contandoli come non trovati, l'AUC scende sotto il caso in gbm (0.413) e a 0.554
  in perio.

## 6. SPATA2 rispetto agli altri metodi

- **AUC sul Golden Standard (universo comune, senza filtro):** SPATA2 è sotto il GLM
  lineare in tutti e sei i dataset.

  | | heme_1000 | stab_inj1 | stab_inj2 | gbm | perio | head_and_neck |
  |---|---|---|---|---|---|---|
  | SPATA2 | 0.931 | 0.676 | 0.698 | 0.722 | 0.793 | 0.708 |
  | GLM lineare | 0.978 | 0.912 | 0.829 | 0.949 | 0.895 | 0.958 |

- **Sovrapposizione dei geni significativi (heme_1000):** dei significativi di
  SPATA2, il GLM lineare ne ritrova il 49.9%, contro il 39.3% atteso a caso.
  Spearman fra i rank: 0.25. Su tutti i dataset lo Spearman fra i rank di SPATA2 e
  del GLM lineare va da 0.04 a 0.50. SPATA2 non è un sottoinsieme più severo del GLM: le due
  statistiche ordinano i geni in modo diverso.
- **Indipendenza dall'espressione:** la correlazione del NES di SPATA2 con il
  controllo sull'espressione è vicina a zero o negativa (da −0.53 a 0.32). Per il GLM
  lineare va da −0.07 a 0.65.
- **Differenze dal tutorial ufficiale:**
  - Il tutorial fa lo screening solo sui geni selezionati come spazialmente variabili
    da SPARK-X (`runSPARKX()`, `getSparkxGenes()`, p < 0.05); la pipeline su tutti i
    geni.
  - Il tutorial usa anche l'adattamento ai modelli di gradiente
    (`@results$model_fits`, MAE < 0.2 come esempio) e ne ricava la direzione
    ("descending", "ascending"); la pipeline usa solo la tabella `significance`
    (FDR e `tot_var`).
  - Il tutorial chiama `identifyTissueOutline()` e `identifySpatialOutliers()` con
    `method = "obs"`; la pipeline con i valori predefiniti.
  - `removeGenesStress()` è facoltativo nel tutorial e non è usato nella pipeline.

---

## 7. Perché l'arricchimento non è un buon benchmark fra metodi

1. **Una classifica banale produce tanti pathway quanto i metodi.** Ordinando i geni
   solo per espressione media, la GSEA dà quanti pathway con FDR < 0.05 dei metodi
   veri. L'ordine casuale ne dà 0.

   | Pathway con FDR < 0.05 (universo comune) | heme_1000 | stab_inj1 | stab_inj2 | gbm | perio | head_and_neck |
   |---|---|---|---|---|---|---|
   | Solo espressione | 750 | 915 | 801 | 1402 | 918 | 1647 |
   | GLM lineare | 827 | 711 | 311 | 1195 | 761 | 1108 |
   | SPATA2 | 133 | 106 | 190 | 348 | 393 | 27 |

   In gbm (`pipeline_v8_noinfl`) l'ordinamento per espressione arriva a 1801.

2. **Anche i pathway di contesto escono dalla sola espressione.** Nei dataset tumorali
   l'ordinamento per espressione trova 107 (gbm) e 209 (head_and_neck) pathway di
   contesto significativi. In head_and_neck sono più dei 190 del GLM.

3. **Il numero di pathway segue la permissività del test, non la qualità
   dell'ordine.** Quando ai GLM si applica un p-value corretto per l'autocorrelazione
   (nullo spettrale), i pathway significativi in heme_1000 passano da 827 a 439, e da
   311 a 48 in stab_inj2. L'AUC sul Golden Standard invece cambia di poco (0.978 →
   0.963 in heme_1000).

4. **Il conteggio dipende da come un metodo filtra i geni.** Con `rm_zero_infl` i
   pathway significativi di SPATA2 in gbm passano da 118 a 205 e la frazione di
   contesto da 0.305 a 0.180, mentre l'AUC di contesto cambia appena (0.662 → 0.634).

5. **La frazione di contesto è calcolata su numeri piccoli per SPATA2.** In
   head_and_neck sono 6 pathway di contesto su 27. Piccole variazioni spostano molto
   la frazione.

6. **L'AUC di contesto è instabile fra dataset e fra metodi.**
   - Per il GLM lineare va da 0.447 (heme_1000) a 0.767 (stab_inj1).
   - La classifica dei metodi cambia da un dataset all'altro: la frazione di contesto
     di SPATA2 è più alta di quella del GLM in 5 dataset su 6 ma non in perio (circa
     6% contro 9%).
   - Il solo ordinamento per espressione ha AUC di contesto da 0.21 a 0.59 secondo il
     dataset.

7. **Il NES segue l'espressione per alcuni metodi.** La correlazione fra il NES del
   GLM lineare e quello del controllo sull'espressione è 0.62 in gbm e 0.55 in
   head_and_neck. Una parte dell'arricchimento del GLM in quei dataset è quindi
   arricchimento per geni molto espressi.

8. **La classifica non ha direzione.** La GSEA è fatta con `scoreType = "pos"` su
   liste ordinate per significatività: i geni che salgono e quelli che scendono con
   la distanza stanno insieme. I pathway in cima riflettono la composizione del
   tessuto più che il processo biologico legato all'annotazione.

9. **I "significativi" non sono calibrati per nessuno dei due metodi principali.**
   - Con un centro finto, scelto lontano dalla lesione, il GLM parametrico trova
     comunque da 1500 a 5000 geni BH-significativi (6 dataset).
   - In heme_1000 un centro qualunque dà circa 3700 geni significativi, contro 6023
     del centro vero.
   - Il nullo di SPATA2 non ha struttura spaziale (§1).

   Contare geni o pathway "significativi" confronta quindi due test entrambi troppo
   ottimistici, ciascuno a modo suo.

10. **Il Golden Standard misura l'ordine e regge.** L'AUC sul Golden Standard usa solo
    l'ordine dei geni, non dipende dalla calibrazione dei p-value e dà un risultato
    coerente in tutti e sei i dataset. Va però tenuto presente che in alcuni dataset il
    Golden Standard è legato all'espressione: l'AUC del solo ordinamento per
    espressione è 0.81 in gbm, 0.87 in head_and_neck, 0.71 in perio, 0.67 in heme_1000,
    e sotto 0.5 nei due stab_inj. In gbm il Golden Standard ha 7280 geni, circa il 35%
    di quelli testati.
