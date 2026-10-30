# Report di sessione: arricchimenti, calibrazione dei test, SPATA2

Report da passare a un'altra conversazione. Riassume cosa è stato fatto sulla
pipeline `pipeline_v7`, cosa è stato scoperto, quali decisioni ha preso
l'autore e cosa resta aperto. Dove un'informazione è una deduzione non
verificata, è indicato.

## 0. Contesto

- **Progetto:** benchmark di metodi per trovare geni con gradiente spaziale in
  trascrittomica spaziale (Visium), rispetto a un centro (sito di lesione o
  iniezione).
- **Metodi confrontati:**
  - GLM radiali (edgeR, lineare e gaussiano sulla distanza dal centro);
  - GLM su centro inferito e su pseudotime (Gaston, SpaceFlow, Slingshot, Monocle);
  - Hybrid_GLM;
  - riferimenti center-free: Moran's I, Laplacian Score, SpatialDE (da CSV);
  - tool esterno SPATA2 (`spatialAnnotationScreening`).
- **Valutazione:**
  - Golden Standard (geni differenziali dell'articolo di ogni dataset): boxplot
    dei rank, ROC/AUC, AUC con intervallo bootstrap;
  - GSEA per metodo.
- **Dataset (6):** heme_1000 (topo, iniezione di eme), stab_inj1 e stab_inj2
  (lesione cerebrale), gbm (glioblastoma, topo), perio (parodontite), head_and_neck
  (tumore). Contesto di arricchimento: "inflammation" per heme, stab, perio;
  "tumor" per gbm e head_and_neck.
- **Codice:** repo GitHub `9nt0n10/RGdist`, branch `claude/nifty-ride-xypva8`
  (nessuna pull request). Il file principale è `pipeline_v7.ipynb` (notebook R).
  Ogni dataset ha un proprio notebook con lo stesso codice comune e un blocco
  dataset diverso: **le modifiche vanno copiate negli altri notebook**.
- **Regola di ordinamento dei geni (decisa dall'autore):** FDR crescente, poi la
  statistica del metodo come spareggio:
  - GLM: `log_p`;
  - Moran: Z;
  - SpatialDE: LLR;
  - SPATA2: `tot_var`.

## 1. Modifiche fatte a `pipeline_v7.ipynb` (stato attuale del branch)

### 1.1 GSEA (sezione 16) semplificata
- Una GSEA per metodo (`clusterProfiler::GSEA`) sulla lista ordinata come sopra.
- Punteggio = posizione (N…1) invece della vecchia scala 1→0. Equivalente per la
  GSEA, ma più semplice.
- `exponent = 0`: KS non pesata, conta solo l'ordine dei geni.
- `scoreType = "pos"`: la classifica non ha segno.
- Nessuna soglia sui pathway: si leggono NES e AUC contesto contro altri sul NES
  (`simple_auc`, intervallo bootstrap, test di permutazione). `p_adjust` resta
  solo come informazione nei CSV.
- Controlli negativi nella stessa GSEA: `Expression_control` (geni ordinati per
  espressione media) e `Random_control`.
- Nuovi grafici: `04_enrichment_context_auc.png` (AUC di contesto con IC),
  `04b` (NES dei 40 pathway di contesto più alti, senza filtro), `04c` (boxplot
  dei NES contesto/altri).
- Salva `enrichment_context_methods.csv` (tabella `enr_tab`).

### 1.2 Sezioni rimosse
- **14d**, delta-AUC appaiato contro Moran's I, e il parametro `AUC_REF_METHOD`.
- **17**, GSEA sui geni comuni a tutti i metodi.
- Sezioni 19 (riepilogo) e 20 (file di verifica) adattate: usano `simple_auc` e
  non hanno più delta-AUC né `enr_delta`.

### 1.3 Nuova regola sui geni: il filtro vale solo per i GLM
Decisione dell'autore: nessun universo comune, ogni metodo con la propria lista
intera, sia per boxplot e AUC sia per la GSEA.
- **Sezione 4:** il filtro "espresso in ≥ `MIN_SPOTS_PER_GENE` (3) spot" resta su
  `spe` e vale solo per i GLM (lineare, gaussiano, inferiti, pseudotime, Hybrid,
  scansione 9c, diagnostiche 10e e 10f). Prima del filtro si salva
  `logcounts_all`: tutti i geni con almeno un conteggio, dopo QC degli spot e
  rimozione dei mitocondriali.
- **10b (SpatialDE):** lista intera del CSV (`sde_full`), ordinata per qval e poi
  LLR. Niente taglio sull'universo, niente coda di non testati.
- **10d (Moran, Laplacian):** calcolati su `logcounts_all`.
- **11b (SPATA2):** lista intera dell'output (`spata2_full`), ordinata per FDR e poi
  `tot_var`, senza coda. Nel riepilogo di SPATA2 ci sono `geni_in_lista` e
  `geni_tolti_dal_filtro_glm`.
- **15:** l'universo per scegliere i pathway e per i controlli è l'unione delle
  liste con almeno un conteggio. `T2_full` = pathway scelti con tutti i loro
  geni, così ogni GSEA li interseca con la propria lista.
- **16:** stampa per metodo `geni` e `tolti_dal_filtro_glm`.
- Tolti i parametri della coda (`SDE_ORIG_NON_TESTATI`, seed della coda).
- **Va rieseguito dalla sezione 4:** i checkpoint vecchi non hanno `logcounts_all`.

Esempio di liste risultanti in perio: GLM 15163 geni; Moran e Laplacian 16790
(1627 espressi in < 3 spot); SpatialDE 15164; SPATA2 8999 con
`rm_zero_infl = TRUE` (vedi §5).

### 1.4 Cosa NON è ancora stato implementato (proposto, in attesa di conferma)
- AUC sul Golden Standard con i geni del GS **non testati da un metodo contati
  come mancati** (messi in fondo). Con liste diverse è l'unica AUC equa.
- Boxplot con **rank relativo** (posizione / lunghezza della lista), GS non
  testati al 100%. Il rank grezzo inganna con liste di lunghezza diversa.
- `METHOD_ORDER` (sezione 12) basato sul rank relativo invece che grezzo.
- Riepilogo con i geni del GS presenti in ogni lista.
- `SPATA2_RM_ZERO_INFL <- FALSE` fisso. Il commento attuale dice "default del
  tool", ma probabilmente il default è FALSE (vedi §5).

## 2. Script autonomi creati (nel repo, fuori dalla pipeline)

| Script | Cosa fa |
|---|---|
| `diagnostica_dataset.R` | legge i CSV di più dataset: % significativi, recupero e fold del GS, GSEA (AUC di contesto, correlazione col controllo sull'espressione), nullo 10f, AUC, calibrazione 10e, confronto SPATA2–altri (sovrapposizione, Spearman, saturazione p), geni significativi di SPATA2 fuori universo |
| `prova_p_empirico*.R` | prova del p empirico dai "centri lontani" (superata, vedi §4) |
| `prova_nulli.R` | file unico: due nulli spaziali per GLM_linear (centri poco correlati, randomizzazione spettrale, più versione sigma e 10000 campi), controllo con falsi centri, GS e GSEA, 6 dataset in parallelo, output in una cartella |
| `confronto_zero_infl.R` | confronta due run della pipeline (SPATA2 con `rm_zero_infl` TRUE e FALSE): liste, significativi, GS (recuperati, AUC propria e con mancanti in fondo, rank relativo), GSEA, dettaglio SPATA2. Controllo di coerenza: cambia solo SPATA2 |

Nessuno script è stato eseguito nell'ambiente di sviluppo (senza R): li ha
eseguiti l'autore, e i numeri sotto vengono dai suoi output.

## 3. Risultati principali

### 3.1 Il numero di geni/pathway significativi non è un buon benchmark
- Il solo ordinamento per espressione dà 750–1800 pathway con FDR < 0.05 per
  dataset, quanti il GLM. Il numero di pathway segue il numero di geni
  significativi e l'espressione.
- L'AUC di contesto del GLM varia molto fra dataset (circa 0.45–0.77).
- La classifica non ha direzione.
- Conclusione condivisa con l'autore: arricchimento come analisi descrittiva,
  non per decidere quale metodo è migliore.

### 3.2 I p-value del GLM sono troppo ottimistici (autocorrelazione spaziale)
- **Significativi:** il GLM lineare parametrico chiama significativo dal 7% al 69%
  dei geni secondo il dataset. Moran (~42%) e SpatialDE (~35%) sono simili in
  heme_1000.
- **Nullo spaziale 10f (heme_1000):**
  - centro vero: 6023 geni significativi (score test);
  - centri lontani: mediana 3707, p empirico 0.0009 (lineare), 0.018 (gaussiano).
  - Il centro conta, ma circa il 60% dei significativi esce con qualunque centro.
- **Calibrazione 10e:** la profondità spiega poco (frazione con p < 0.05 sotto il
  nullo 0.06–0.17, contro 0.45–0.55 osservata).
- **Falsi centri (`prova_nulli.R`, 6 dataset):** con un centro finto il GLM
  parametrico trova 1500–5000 geni BH-significativi, e il 20–36% dei geni ha
  p < 0.05.

### 3.3 Correzioni provate per il p del GLM
- **Centri lontani** (oltre la mediana della distanza): difettosa. Le distanze
  sono correlate con quella vera; la versione con segno era più permissiva del
  parametrico in stab_inj1, stab_inj2 e head_and_neck.
- **Centri con |r| < 0.3:** non calibrata. Con i falsi centri dà ancora
  1400–3600 geni BH-significativi (0 solo in perio).
- **Randomizzazione spettrale:** campi casuali con lo stesso spettro di
  autocorrelazione del raggio, sugli autovettori della matrice dei 6 vicini.
  - Moran's I dei campi circa 1.0, uguale al raggio.
  - **È l'unica calibrata:** falsi centri con 0 geni BH-significativi in 5 dataset
    su 6 (gbm 204, contro 5081 del parametrico).
  - Costo: AUC sul GS −0.01/−0.08.
  - stab_inj1 e stab_inj2: 0 significativi anche col centro vero (risoluzione con
    1000 campi).
  - Preparate ma non lanciate: versione sigma, p = 2Φ(−|Z|/σ_gene), e versione
    con 10000 campi.
- **Decisione dell'autore:** tenere il GLM parametrico per l'ordine dei geni (le
  ROC peggioravano con le correzioni). Il nullo spettrale serve solo se si
  vogliono presentare conteggi di geni significativi.

### 3.4 L'ordine dei geni del GLM è solido, ma in parte legato all'espressione
- AUC sul GS dei GLM con centro: 0.83–0.98, la migliore o fra le migliori in ogni
  dataset. SPATA2 è sotto il GLM in tutti e sei.
- AUC sul GS ordinando solo per espressione media: gbm 0.81, head_and_neck 0.87,
  perio 0.71, heme 0.67, stab_inj 0.33–0.40. In alcuni dataset il GS è legato
  all'espressione.
- In gbm il GS ha 7280 geni (circa 35% dei testati): sopra il 30% la pipeline
  avverte che boxplot e ROC distinguono poco.

### 3.5 SPATA2: cosa fa davvero (dal paper Nat Commun 2024 e dal sito)
- **p-value:** frazione di 10000 profili simulati da "valori casuali assegnati a
  tutti gli spot" con `tot_var` ≤ quella osservata. Il nullo non ha struttura
  spaziale: stesso problema del GLM. Nel paper, un esempio sul glioblastoma dà
  11.478 geni significativi.
- **Saturazione:** in gbm (senza filtro) 3450 geni hanno p = 0. Fra loro decide solo
  `tot_var`. Probabilmente il p è funzione della sola `tot_var` (nullo uguale per
  tutti i geni): **da verificare** con Spearman fra `p_value` e `tot_var`.
- **Artefatto sui geni sparsi:** fra le posizioni 1500 e 2500 di SPATA2 in gbm,
  891/1000 geni sono espressi in < 3 spot, con `tot_var` = 2.000 e p = 0. È il
  "buco" senza geni del GS nel boxplot. Spiegazione dedotta: SPATA2 riscala il
  profilo fra 0 e 1, e un gene acceso in uno spot dà una punta 0→1→0 con
  variazione totale 2. In heme_1000, senza filtro, SPATA2 chiamava significativo
  il 96% (1681/1747) dei geni sparsi.
- **Preprocessing nella pipeline (11b) rispetto al tutorial ufficiale:**
  - uguali: `initiateSpataObjectVisium`, `removeGenesZeroCounts`,
    `removeObsZeroCounts`, `normalizeCounts` (LogNormalize), `core = FALSE`,
    distanza predefinita (bordo del tessuto);
  - diversi: `identifyTissueOutline` e `identifySpatialOutliers` sono chiamati
    con i default, mentre il tutorial usa `method = "obs"`, `eps = "125um"`,
    `minPts = 3`;
  - manca `removeGenesStress()` (facoltativo nel tutorial);
  - **differenza principale:** il tutorial fa lo screening solo sui geni
    selezionati da SPARK-X (`runSPARKX`, `getSparkxGenes`, p < 0.05), la
    pipeline su tutti i geni. Per il benchmark serve una classifica completa,
    quindi è uno scostamento da dichiarare;
  - il tutorial usa anche `@results$model_fits` (MAE < 0.2, come esempio) e i
    modelli "descending"/"ascending" per la **direzione**: SPATA2 può dare un
    segno, che la pipeline non usa.

### 3.6 SPATA2 e `rm_zero_infl`
- **Regola dal paper:** si tolgono gli outlier oltre 1.5 IQR; se restano solo zeri,
  il gene è "zero-inflated" ed escluso. In pratica (deduzione) toglie i geni
  rilevati in meno di circa il 25% degli spot. Il paper dice che l'opzione "can be
  activated" ponendola a TRUE, quindi il default sembra FALSE (da verificare con
  `formals(SPATA2::spatialAnnotationScreening)$rm_zero_infl`). Gli autori ne
  consigliano l'uso.
- **perio, con TRUE:**
  - lista di 8999 geni contro 15163;
  - perde 126/353 geni del GS, quasi tutti dell'infiltrato immunitario (*CD3D*,
    *LCK*, *PTPRC*, *CXCR4*, *CCR7*, *CD79A*, immunoglobuline);
  - AUC SPATA2: 0.861 sulla propria lista, 0.554 con i mancanti in fondo;
  - GLM lineare: 0.895 (0.915–0.927 per i GLM con centro sui soli geni di SPATA2).
- **gbm, TRUE contro FALSE (`confronto_zero_infl.R`, run `benchmark_v8` e
  `benchmark_v8_noinfl`):**
  - lista 5455 contro 20503 (−73%);
  - GS nella lista 3549 contro 7280; tolti anche marcatori di proliferazione
    (*Mki67*, *Cdk1*, *Plk1*, *Bub1*);
  - recall GS 30.2% contro 50.7%;
  - AUC propria 0.847 contro 0.633; con i mancanti in fondo 0.413 contro 0.633;
  - ordine dei geni rimasti identico (Spearman 1);
  - GLM lineare 0.949 in entrambe le letture;
  - controllo di coerenza ok (cambia solo SPATA2).
- **Conclusione:** con `rm_zero_infl = TRUE` l'AUC di SPATA2 sale solo perché
  metà del GS non viene valutata. Per il benchmark: FALSE. Un gene non testato
  equivale a "nessuna evidenza", quindi va in fondo; non si può dargli un p-value
  calcolato con un altro test senza smettere di valutare SPATA2.
- **Copertura del GS senza filtro (gbm):** SPATA2 ha 7280 geni del GS, SpatialDE
  7275, GLM, Moran e Laplacian 7267. Nessun metodo li ha tutti: geni non rilevati
  nel campione, filtri propri, simboli non corrispondenti.

### 3.7 Dove SPATA2 va meglio del GLM
- Sull'arricchimento di contesto. In gbm (liste intere):
  - AUC di contesto 0.66 (0.63 con filtro) contro 0.57 del GLM lineare;
  - quota di pathway di contesto fra i significativi 0.31 contro 0.14;
  - correlazione del NES con l'espressione −0.07 contro 0.65.
- Nella versione precedente (universo comune) la quota di contesto era più alta
  per SPATA2 in 5 dataset su 6 (non in perio). Ma SPATA2 trova pochi pathway
  significativi (27–393), quindi la quota è calcolata su numeri piccoli.

## 4. Decisioni dell'autore (da rispettare)
1. Ordinamento per FDR più spareggio del metodo; GSEA solo sull'ordine
   (`exponent = 0`).
2. Nessuna soglia sui pathway nella GSEA.
3. Niente confronto appaiato con Moran (14d) né GSEA sui geni comuni (17).
4. Filtro sui geni solo per i GLM; ogni altro metodo con la sua lista intera;
   nessun universo comune e nessuna coda.
5. GLM con p-value parametrico (le correzioni spaziali peggioravano le ROC).
6. Orientamento: arricchimento poco adatto come benchmark fra metodi; Golden
   Standard come misura principale.

## 5. Punti aperti
- Implementare quanto in §1.4: AUC con GS mancanti in fondo, rank relativo,
  `METHOD_ORDER`, copertura del GS, `rm_zero_infl = FALSE` con commento corretto.
- Verificare che il p di SPATA2 dipenda solo da `tot_var`.
- Allineare il preprocessing di SPATA2 al tutorial (`method = "obs"` ecc.) e
  decidere su `removeGenesStress()` e SPARK-X (con i geni non selezionati da
  trattare come non testati).
- Eventuale prova "falsi centri" per SPATA2 (annotazioni finte) per misurarne i
  falsi positivi.
- Nullo spettrale con sigma o 10000 campi, solo se servono conteggi di
  significativi.
- Ripetere il confronto `rm_zero_infl` sugli altri dataset.
- Eventuale GSEA con segno (GLM lineare dal coefficiente, SPATA2 dai modelli).
- Copiare il codice comune nei notebook degli altri cinque dataset.
- Laplacian_Score appare con 0% significativi: va mostrato come NA (non ha
  p-value).
- Ricalcolo dell'FDR di SPATA2: superato dalla regola delle liste intere (ora
  l'FDR è quello di SPATA2 sul suo output).

## 6. Messaggio per i colleghi (versione concordata)

> Ho lavorato un po' sugli arricchimenti di SPATA2 e ho scoperto alcune cose su
> cui avrei bisogno di un consiglio. Innanzitutto avevo commesso un errore a
> favore di SPATA2: gli avevo applicato il filtro che uso per i GLM, mentre lui
> ha già i suoi filtri. Corretto questo, SPATA2 è peggiorato ulteriormente. Uno
> dei motivi è che con i geni con molti zeri il tool va un po' in tilt (da qui il
> buco nel boxplot dei rank: i primi migliaia di geni hanno tutti p-value 0 e
> sono ordinati per la tot_var, che uso come spareggio, così i geni espressi in
> pochissimi spot finiscono in alto). Per questo gli autori consigliano di
> attivare un parametro, `rm_zero_infl`. L'ho provato (in gbm e in parodontite) e
> sembra migliorare molto, troppo: i rank dei geni del Golden Standard diventano
> bassissimi. Guardando meglio, è perché il parametro scarta i geni rilevati in
> pochi spot, cioè una gran parte dei geni che gli altri tool testano, e di cui
> quindi non abbiamo nemmeno il p-value. Scarta così anche circa metà dei geni
> del Golden Standard in gbm e un terzo in parodontite. Contando quei geni come
> non trovati, l'AUC di SPATA2 crolla. La frazione di pathway di contesto
> invece resta quasi sempre a favore di SPATA2 (5 dataset su 6), ma credo che
> gli arricchimenti non siano un buon benchmark, anche perché SPATA2 trova pochi
> pathway significativi. Cosa ne pensate? Lascereste perdere gli arricchimenti
> come confronto fra metodi e terreste solo il Golden Standard?
