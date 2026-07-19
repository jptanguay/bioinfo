# L'algorithme de de Bruijn et la recherche de chemins eulériens : de la bioinformatique aux architectures réseau

## Introduction

En 1946, le mathématicien hollandais Nicolaas de Bruijn cherche à résoudre un problème de logique pure : comment construire la plus courte suite de symboles contenant toutes les combinaisons possibles d'une longueur donnée. Pour y parvenir, il imagine un réseau géométrique de mots reliés entre eux. Cet objet mathématique prendra le nom de **graphe de de Bruijn**.

Pendant des décennies, cette découverte reste confinée aux livres de mathématiques discrètes. Tout change au début des années 2000 avec l'explosion du séquençage de l'ADN. Confrontés à des milliards de courts fragments d'ADN à remettre dans l'ordre, les bioinformaticiens ont réalisé que le graphe de de Bruijn était l'outil idéal pour reconstruire informatiquement les génomes.

Ce chapitre explore le fonctionnement de cet algorithme, la manière dont les ordinateurs résolvent ce puzzle à l'aide de parcours de graphes, et comment cette même structure mathématique aide aujourd'hui à concevoir des superordinateurs ou à guider des robots.

## 1. Le défi du séquençage et de l'assemblage d'ADN

### 1.1. Le problème biologique : les fragments de lecture (*reads*)

Le génome d'un organisme est une immense chaîne textuelle composée de quatre bases chimiques représentées par les lettres A, T, C et G. Chez l'humain, cette chaîne s'étend sur environ 3 milliards de caractères.

Les technologies de séquençage actuelles (le "séquençage de nouvelle génération" ou NGS) font face à une contrainte physique majeure : elles sont incapables de lire un chromosome entier d'un seul coup. Pour contourner cette limite, les biologistes cassent l'ADN d'un échantillon en des millions de copies identiques, que la machine lit simultanément. Le résultat obtenu est un fichier brut contenant des millions de courtes séquences textuelles appelées **reads**, mesurant généralement entre 100 et 150 lettres.

Le travail du bioinformaticien consiste à faire un **assemblage de novo**, c'est-à-dire reconstruire le texte d'origine uniquement à partir de ces *reads*, sans modèle pour s'aider.

```
Texte d'origine inconnu :  [?????????????????????????????????]

Reads obtenus :            [ATGCTA]
                              [GCTAAA]
                                 [TAAAGC]
                                    [AAGCTT]
```

### 1.2. La notion de $k$-mer

Pour rendre les données exploitables par un ordinateur, la première étape consiste à découper chaque *read* en sous-chaînes de longueur fixe appelée $k$. Un fragment de longueur $k$ est appelé un **$k$-mer**.

L'extraction des $k$-mers se fait par une technique de fenêtre glissante : on avance d'une lettre à la fois le long du *read*.

> **Exemple concret :**
> 
> Soit un *read* égal à `ATGCTA` et une taille $k = 4$.
> 
> - Fenêtre 1 : `ATGC`
> 
> - Fenêtre 2 : `TGCT`
> 
> - Fenêtre 3 : `GCTA`
> 
> Le fragment `ATGCTA` génère donc **trois** 4-mers.

Le choix de la variable $k$ est crucial. Si $k$ est trop petit, le graphe sera très compact mais saturé de fausses liaisons. Si $k$ est trop grand, la moindre erreur de lecture éliminera des données valides.

### 1.3. Les limites des anciennes approches (OLC)

Aux débuts de la bioinformatique, pour le projet du génome humain, les logiciels utilisaient la méthode **Overlap-Layout-Consensus (OLC)**. Cette méthode consistait à comparer tous les *reads* deux à deux pour trouver leurs zones de chevauchement (*overlap*), puis à dessiner un plan d'agencement (*layout*).

Le problème de la méthode OLC est sa complexité informatique. Comparer chaque fragment avec tous les autres requiert un nombre d'opérations proportionnel au carré du nombre de fragments, noté mathématiquement $O(N^2)$. Avec l'arrivée des séquenceurs modernes produisant des centaines de millions de *reads*, le temps de calcul requis et la quantité de mémoire vive nécessaire deviennent si gigantesques que la méthode OLC est inapplicable en pratique. Pour surmonter ces contraintes de complexité, l'approche par graphe de de Bruijn a été introduite dans les algorithmes d'assemblage, notamment par Pavel Pevzner en 2001, modifiant ainsi le paradigme du traitement des données de séquençage.

## 2. Structure et construction du graphe de de Bruijn

### 2.1. Définition mathématique du graphe

Un graphe de de Bruijn est un réseau orienté construit à partir d'un alphabet de symboles (dans notre cas : A, T, C, G). Sa principale particularité réside dans l'identité de ses éléments :

- Les **sommets (ou nœuds)** du graphe représentent des séquences de longueur $k-1$.

- Les **arcs** du graphe représentent les **$k$-mers** réels trouvés dans les données de séquençage.

De manière générale, pour une taille de $k$-mer donnée, les sommets du graphe représentent des séquences de longueur $k-1$. Un arc orienté relie un sommet $A$ à un sommet $B$ si et seulement si le suffixe de longueur $k-2$ du sommet $A$ est identique au préfixe de longueur $k-2$ du sommet $B$. L'arc ainsi formé représente le $k$-mer unique obtenu par la fusion de ces deux sommets.

### 2.2. Algorithme de construction pas à pas

La construction du graphe s'effectue par un traitement séquentiel de l'ensemble des $k$-mers préalablement extraits des *reads* (voir section 1.2). Pour chaque $k$-mer :

1. On isole ses $k-1$ premières lettres : c'est le **préfixe**.

2. On isole ses $k-1$ dernières lettres : c'est le **suffixe**.

3. On crée un nœud pour le préfixe et un nœud pour le suffixe (s'ils n'existent pas déjà dans le graphe).

4. On trace un arc orienté qui part du nœud préfixe et qui pointe vers le nœud suffixe. L'étiquette de cette flèche est le $k$-mer lui-même.

```
Pour le 4-mer "ATGC" :
Nœud Préfixe (3-mer) : [ATG]
Nœud Suffixe (3-mer) : [TGC]

Graphe généré :        [ATG] -------> [TGC]
```

Si un même $k$-mer est détecté plusieurs fois dans les données de séquençage, aucun nouvel arc n'est ajouté au graphe. L'algorithme se contente d'incrémenter le **poids** (ou la multiplicité) de l'arc existant. Par conséquent, la taille du graphe de de Bruijn dépend uniquement de la complexité intrinsèque du génome séquencé, et non de la quantité brute de données générées par le séquenceur.

### 2.3. Modélisation : le chemin eulérien

Une fois que tous les *reads* ont été convertis, le graphe contient l'ensemble de l'information génétique sous forme dispersée. Reconstruire la séquence de l'ADN d'origine consiste à trouver un parcours qui passe par tous les arcs du graphe (puisque chaque arc est un morceau d'ADN valide qu'il faut intégrer).

En théorie de graphe, un parcours qui visite chaque arc d'un graphe exactement une et une seule fois s'appelle un **chemin eulérien**. Si le parcours revient à son point de départ, on parle de **cycle eulérien**. L'assemblage de l'ADN est donc la résolution informatique d'un problème de chemin eulérien.

### 2.4. Les bruits de la réalité biologique

Dans un monde théorique parfait, le graphe de de Bruijn forme une ligne continue facile à lire. En réalité, deux phénomènes biologiques perturbent la structure du graphe :

- **Les erreurs de séquençage :** Les technologies de séquençage introduisent inévitablement du bruit dans les données, principalement sous forme d'erreurs de substitution (par exemple, la détection erronée d'une guanine à la place d'une adénine). Ces artefacts génèrent de faux $k$-mers de faible fréquence qui se traduisent dans le graphe par des structures en cul-de-sac (*dead-ends*) ou par des chemins alternatifs parallèles de courte longueur, appelés "bulles". Les logiciels d'assemblage éliminent ces structures parasites via des algorithmes d'élagage, en supprimant les chemins dont le poids (la couverture) est anormalement bas par rapport à la moyenne du graphe.

- **Les répétitions génomiques :** Les génomes contiennent de nombreuses séquences répétitives dont la longueur dépasse celle du $k$-mer choisi (par exemple, des motifs tels que `ATATAT...` peuvent se répéter plusieurs fois de manière contiguë). 
  
  Dans le graphe, ces régions identiques convergent vers un sous-graphe commun, créant des intersections complexes. Lors de la recherche du chemin eulérien, l'algorithme est confronté à une ambiguïté de routage : il ne dispose pas d'assez d'informations locales pour déterminer quel arc sortant correspond logiquement à quel arc entrant. 
  
  Pour résoudre ces bifurcations, des méthodes avancées de **désambiguïsation** (ou *scaffolding*) sont employées. Ces techniques s'appuient sur des données de contexte à plus longue portée, telles que le séquençage par paires (*paired-end*), les technologies de lectures longues (*long reads*), ou encore la cartographie optique du génome. Ces méthodes de résolution d'ambiguïtés dépassent toutefois le cadre de ce chapitre et ne seront pas détaillées ici.

## 3. Algorithmes de recherche du chemin eulérien

### **3.1. Conditions d'existence du chemin**

Avant de lancer les calculs, l'algorithme doit vérifier qu'un chemin eulérien existe bien dans le graphe. Pour un graphe orienté d'un seul tenant (connexe), les règles de validation sont simples et reposent sur le décompte des arcs reliés à chaque sommet :

- **Tous les sommets** doivent avoir un **degré entrant égal au degré sortant**,  **sauf éventuellement deux sommets** :
  
  - **Un sommet** a un **degré sortant = degré entrant + 1** (c’est le sommet de départ du chemin).
  - **Un sommet** a un **degré entrant = degré sortant + 1** (c’est le sommet d’arrivée du chemin).

- Si **tous les sommets** ont un degré entrant égal au degré sortant, alors le chemin eulérien est un **cycle eulérien** (il commence et finit au même sommet).

Si le graphe respecte ces conditions, l'existence du chemin est garantie et l'algorithme peut commencer la reconstruction.

### 3.2. L'algorithme de Hierholzer

L'algorithme de Hierholzer est la méthode privilégiée en programmation bioinformatique. Sa force réside dans sa rapidité : sa complexité temporelle est linéaire, soit $O(M)$ où $M$ est le nombre total d'arc. L'algorithme ne visite chaque arc qu'une seule fois.

#### Principe de fonctionnement

L'algorithme repose sur le concept de fusion de sous-cycles. Il utilise une structure de pile informatique pour mémoriser le cheminement actuel.

```
Fonction Trouver_Chemin_Hierholzer(Graphe G, Noeud depart):
    Pile_Travail = [depart]
    Chemin_Final = []

    Tant que Pile_Travail n'est pas vide:
        Noeud_Courant = Regarder_Sommet(Pile_Travail)

        Si G contient des arcs sortants depuis Noeud_Courant:
            Prochain_Noeud = Choisir_Voisin_Et_Supprimer_Arc(G, Noeud_Courant)
            Ajouter Prochain_Noeud à Pile_Travail
        Sinon:
            Noeud_Bloque = Retirer_Sommet(Pile_Travail)
            Ajouter Noeud_Bloque à Chemin_Final

    Inverser_Ordre(Chemin_Final)
    Retourner Chemin_Final
```

L'algorithme avance à l'aveugle dans le graphe. Dès qu'il se retrouve coincé dans une impasse, il sait qu'il a atteint la fin d'un sous-cycle. Il enregistre ce nœud dans la liste finale, puis "recule" d'un pas pour explorer les autres flèches restantes. Les morceaux de cycles se recollent ainsi automatiquement dans le bon ordre.

### 3.3. L'algorithme de Fleury

L'algorithme de Fleury est une approche alternative, plus intuitive pour un être humain mais inefficace pour un ordinateur. Sa complexité est de $O(M^2)$, ce qui signifie que doubler la taille du graphe multiplie le temps de calcul par quatre.

#### Principe de fonctionnement

Fleury progresse de manière frontale en effaçant les acrs derrière lui au fur et à mesure qu'il avance. Sa seule règle est de ne jamais traverser un **pont** (un arc dont la suppression couperait le graphe restant en deux parties isolées), sauf s'il n'y a absolument aucun autre chemin disponible.

```
Fonction Trouver_Chemin_Fleury(Graphe G, Noeud depart):
    Noeud_Courant = depart
    Chemin = [depart]

    Tant que G contient des arcs:
        Arcs_Dispo = Arcs sortant de Noeud_Courant
        Choix = Arcs_Dispo[0]

        Pour chaque F dans Arcs_Dispo:
            Si NON Verifier_Si_Pont(G, F):
                Choix = F
                Quitter Boucle Pour

        Noeud_Courant = Destination_De(Choix)
        Ajouter Noeud_Courant à Chemin
        Supprimer_Arc_Du_Graphe(G, Choix)

    Retourner Chemin
```

#### Pourquoi Fleury échoue en pratique

La fonction `Verifier_Si_Pont` oblige l'ordinateur à simuler la suppression de l'arc et à recalculer l'accessibilité de l'intégralité du graphe restant (via un parcours en largeur ou en profondeur). Répéter cette vérification globale à chaque pas de progression demande une puissance de calcul colossale, rendant Fleury inutilisable sur des données de taille génomique.

## 4. Au-delà de la biologie : Autres applications des graphes de de Bruijn

Au-delà de l'assemblage de génomes ADN *de novo*, les propriétés combinatoires des graphes de de Bruijn les rendent incontournables dans plusieurs autres domaines de la bioinformatique :

- **Transcriptomique (Séquençage d'ARN / RNA-Seq) :** Utilisation des graphes pour assembler le transcriptome afin d'identifier les différents transcrits et caractériser le phénomène d'épissage alternatif (*alternative splicing*).

- **Métagénomique :** Assemblage et analyse de mélanges complexes d'ADN issus de communautés microbiennes environnementales ou cliniques (comme le microbiote intestinal), sans nécessiter l'isolement préalable des espèces.

- **Détection de variants (Variant Calling) :** Identification directe de mutations, d'insertions ou de délétions en comparant la topologie du graphe de de Bruijn des données séquencées à celle d'un génome de référence.

- **Compression de données de séquençage :** Optimisation du stockage des fichiers de *reads* massifs en exploitant la structure compacte du graphe pour éliminer la redondance des séquences.

Voici la section documentaire pour clore ce chapitre. Elle regroupe les publications scientifiques fondatrices ainsi que les ressources en ligne de référence pour approfondir l'étude des graphes de de Bruijn et leurs applications.

## 5. Ressources et littérature scientifique de référence

### 5.1. Publications académiques fondatrices

- **L'article mathématique originel :**
  
  - De Bruijn, N. G. (1946). *A combinatorial problem*. Koninklijke Nederlandse Akademie van Wetenschappen, 49, 758-764.
  
  - *Note :* C'est le texte historique où le mathématicien introduit la structure de graphe pour résoudre le problème des suites de symboles.

- **L'introduction des graphes de de Bruijn en bioinformatique :**
  
  - Pevzner, P. A., Tang, H., & Waterman, M. S. (2001). *An Eulerian path approach to DNA fragment assembly*. Proceedings of the National Academy of Sciences (PNAS), 98(17), 9748-9753.
  
  - *Note :* Cet article capital marque le passage de la méthode OLC aux graphes de de Bruijn pour le séquençage de nouvelle génération (NGS).

- **Les algorithmes des assembleurs modernes :**
  
  - Zerbino, D. R., & Birney, E. (2008). *Velvet: Algorithms for de novo short read assembly using de Bruijn graphs*. Genome Research, 18(5), 821-829.
  
  - Bankevich, A., et al. (2012). *SPAdes: A new genome assembly algorithm and its applications to single-cell sequencing*. Journal of Computational Biology, 19(5), 455-477.
  
  - *Note :* Ces deux publications décrivent concrètement comment les logiciels intègrent la théorie des graphes pour nettoyer les erreurs biologiques (bulles et impasses).

### 5.2. Liens web et ressources pédagogiques fiables

- **Le portail de référence NCBI (National Center for Biotechnology Information) :**
  
  - [ncbi.nlm.nih.gov](https://www.ncbi.nlm.nih.gov/) — Le site officiel de la plus grande base de données de bioinformatique au monde. On y trouve la section *Bookshelf* qui héberge des manuels complets en accès libre sur le traitement des données de séquençage et l'assemblage de génomes.

- **Cours en ligne de référence (Rosalind Project) :**
  
  - [rosalind.info](http://rosalind.info/) — Une plateforme éducative spécialisée et gratuite conçue pour apprendre la bioinformatique par la pratique. La section *Bioinformatics Textbooks Track* contient des fiches explicatives interactives remarquables sur la construction des graphes de de Bruijn et le code de l'algorithme de Hierholzer.

- **Ressources algorithmiques de l'Interstices (INRIA) :**
  
  - [interstices.info](https://interstices.info/) — Ce portail de médiation scientifique, soutenu par l'Institut national de recherche en informatique et en automatique, propose des articles de vulgarisation de haut niveau écrits par des chercheurs, expliquant de manière visuelle la théorie des graphes et la notion de complexité algorithmique.

## Lexique

- **$k$-mer :** Sous-séquence textuelle continue de longueur fixe $k$, extraite par découpage glissant le long d'une séquence de caractères plus longue.

- **Arcs (ou arêtes) :** Liens orientés (flèches) reliant les nœuds d'un graphe. Dans le graphe de de Bruijn appliqué à l'ADN, un arc correspond à un $k$-mer.

- **Assemblage de novo :** Méthode bioinformatique consistant à reconstruire l'intégralité de la séquence d'un génome inconnu à partir de fragments courts, sans s'aider d'une séquence de référence préexistante.

- **Chemin eulérien :** Parcours continu au sein d'un graphe qui visite chaque arc (ou arête pour un graphe non orienté) exactement une seule fois.

- **Complexité algorithmique :** Mesure mathématique (souvent notée sous la forme $O$ majuscule, "Big O") évaluant l'évolution du temps de calcul ou de l'espace mémoire requis par un programme informatique lorsque la quantité de données augmente.

- **Graphe orienté :** Structure mathématique composée de sommets reliés par des flèches possédant un sens de circulation obligatoire (de la source vers la destination).

- **Pont (ou arête de coupure) :** Arête d'un graphe dont la suppression a pour conséquence mathématique d'augmenter le nombre de composants déconnectés du graphe (séparation en deux morceaux isolés).

- **Reads (ou fragments de lecture) :** Courtes séquences de caractères (généralement composées de 100 à 150 bases A, T, C, G) générées directement par les machines de séquençage lors de l'analyse de l'ADN.

- **Sommets (ou nœuds) :** Points d'intersection ou de jonction au sein d'un graphe. Dans le contexte de de Bruijn, ils représentent les préfixes ou suffixes de taille $k-1$.

- **Suite de de Bruijn :** Séquence cyclique de symboles issue d'un alphabet donné dans laquelle chaque combinaison possible de longueur $n$ apparaît exactement une fois et une seule.
