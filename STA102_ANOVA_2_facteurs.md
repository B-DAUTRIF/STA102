---
title: "Analyse de variance (ANOVA) à 2 facteurs"
output: html_notebook
---

Ce bloc-note reprend le contenu du cours STA102 du Conservatoire National des Arts et Métier (CNAM) sur les modéles linéaires. Son but est de décrire les étapes d'une analyse de la variance (ANOVA) à deux facteur, en rassemblant à la fois le matériel théorique et pratique pour travailler efficacement sur des données.

Ce script est destiné à un usage personnel. J'essaie au maximum qu'il puisse être utilisé tel quel pour effectuer les calculs et la production des différents résultats (matrices, vecteurs paramètres et de variance, testes statistiques) : a condition pour cela d'effectuer les changement nécessaires pour l'import de nouveaux jeux de données.

Un éventuel lecteur intéressé par ce document est invité à laisser suggestions et remarques (par le biais de GitHub).

### Présentation des données :

On dispose de $n$ observations sur deux variables explicatives ($X_{1} \ ; \ X{2}$) et une variable **quantitative** qu'on souhaite expliquer ($Y$) :

- On note $i$ le niveau du facteur $X_{1} \ : \ \{x_{1i},...,x_{pi}\}$        
- On note $j$ le niveau du facteur $X_{2} \ : \ \{x_{1i},...,x_{qj}\}$        
- $\forall (i \ ; \ j) \space : \space Y=y_{ij}$

Par convention, on appel chaque couple $(i;j)$ des niveaux de facteurs des **traitements**.

On a donc : 

- $n_{ij}$ : le nombre de répétitions par traitement        
- $k$ : l'indice de répétition du couple $(i;j)$

### Table de contingence : 

On calcul d'abord les effectifs marginaux pour chaque traitement :

$$n_{i.}=\sum_{j=1}^{q}y_{ij}$$

$$n_{.j}=\sum_{i=1}^{p}y_{ij}$$

Puis l'effectif total :

$$n=\sum_{i=1}^{p}\sum_{j=1}^{q}n_{ij}=\sum_{j=1}^{q}n_{i.}=\sum_{i=1}^{p}n_{.j}$$

Puis on produit le tableau des fréquence pour chaque traitement (appelé table de contingence) :

|       | $X_1=1$ | $\cdots$ | $X_1=j$ | $\cdots$ | $X_1=q$ | **Tot.** |
|:-----:|:-------:|:--------:|:-------:|:--------:|:-------:|:--------:|
| $X_2=1$ | $n_{11}$ | $\cdots$ | $n_{1j}$ | $\cdots$ | $n_{1q}$ | $n_{1\cdot}$ |
| $\vdots$ | $\vdots$ |          | $\vdots$ |          | $\vdots$ | $\vdots$ |
| $X_2=i$ | $n_{i1}$ | $\cdots$ | $n_{ij}$ | $\cdots$ | $n_{iq}$ | $n_{i\cdot}$ |
| $\vdots$ | $\vdots$ |          | $\vdots$ |          | $\vdots$ | $\vdots$ |
| $X_2=p$ | $n_{p1}$ | $\cdots$ | $n_{pj}$ | $\cdots$ | $n_{pq}$ | $n_{p\cdot}$ |
| **Tot.** | $n_{\cdot 1}$ | $\cdots$ | $n_{\cdot j}$ | $\cdots$ | $n_{\cdot q}$ | $n$ |

### Moyennes : 

**Rappel :** On note $y_{ijk}$, l'individu $k$ ayant reçu le traitement $ij$. 

On dispose donc de : 

- $Y_{ij.}$ : la moyenne de $Y$ pour le traitement $(i;j)$ :

$$Y_{ij.}=\frac{1}{n_{ij}}\sum_{k=1}^{n_{ij}}y_{ijk}$$

- $Y_{.j.}$ : la moyenne de $Y$ pour le niveau $j$ du second facteur :

$$Y_{.j.}=\frac{1}{n_{i.}}\sum_{i=1}^{p}\sum_{k=1}^{n_{ij}}y_{ijk}$$

- $Y_{i..}$ : la moyenne de $Y$ pour le niveau $i$ du premier facteur :

$$Y_{i..}=\frac{1}{n_{i.}}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}y_{ijk}$$

- $Y_{...}$ : la moyenne générale de $Y$ :

$$Y_{...}=\frac{1}{n_{i.}}\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}y_{ijk}$$

La table de contingence s'organise comme suit :

|       | $X_1=1$ | $\cdots$ | $X_1=j$ | $\cdots$ | $X_1=q$ | **Tot.** |
|:-----:|:-------:|:--------:|:-------:|:--------:|:-------:|:--------:|
| $X_2=1$ | $Y_{11}$ | $\cdots$ | $Y_{1j}$ | $\cdots$ | $Y_{1q}$ | $Y_{1\cdot}$ |
| $\vdots$ | $\vdots$ |          | $\vdots$ |          | $\vdots$ | $\vdots$ |
| $X_2=i$ | $Y_{i1}$ | $\cdots$ | $Y_{ij}$ | $\cdots$ | $Y_{iq}$ | $Y_{i\cdot}$ |
| $\vdots$ | $\vdots$ |          | $\vdots$ |          | $\vdots$ | $\vdots$ |
| $X_2=p$ | $Y_{p1}$ | $\cdots$ | $Y_{pj}$ | $\cdots$ | $Y_{pq}$ | $Y_{p\cdot}$ |
| **Tot.** | $Y_{\cdot 1}$ | $\cdots$ | $Y_{\cdot j}$ | $\cdots$ | $Y_{\cdot q}$ | $Y_{\cdot \cdot}$ |

### Plan d'expérience :

*L'ensemble des $n_{ij}$ définit le plan d'expérience*.

Pour étudier simplement l'effet de chacun des facteurs ($X_{1} ; X_{2}$), le plan d'expérience doit être **orthogonal** :

Définitions : un plan d'expérience est :

- complet si : $\forall \space (i \space ; j) \space : n_{ij} > 0$          
- à répétition si : $\forall \space (i \space ; j) \space : n_{ij} > 1$          
- équilibré si : $\forall \space (i \space ; j) \space : n_{ij} =r > 0$          
- orthogonal si : $\forall \space (i \space ; j) \space : n_{ij}=\frac{n_{i.}\times n_{.j}}{n}$          

Remarque : Tout plan équilibré est orthogonal.

**Attention :** Ne seront abordés ici que les plans d'expérience **équilibrés avec ou sans répétition** ($r\geq 1$).

### Modèle d'analyse de la variance à 2 facteurs :

#### Modèle 1 :

$$Y_{ijk}=\mu_{ij}+\epsilon_{ijk} \space ; \space \epsilon \space iid \sim \mathcal{N}(\mu \space ;\space \sigma^{2})$$

où :

- $\mu_{ij}$ est la moyenne pour chaque traitement.                  
- $\sigma^{2}$ est la variance intra-groupe, commune à tous les traitements.                  

#### Modèle 2 :

$$Y_{ijk}=\mu+\alpha_{i}+\beta_{j}+\gamma_{ij}+\epsilon_{ijk} \space ; \space \epsilon \space iid \sim \mathcal{N}(\mu \space ;\space \sigma^{2})$$

où :

- $\mu$ est un terme constant, représentant la moyenne générale (effet global inconnu) commune à tous les traitements.                  
- $\{\alpha_{i} \space ; \space \beta_{i}\}$ sont les effets de chaque traitement.              
- $\gamma_{i}$ ets un terme d'interraction entre les traitements. C'est lui qui permet de décrire au delà d'une somme les différents effets.                       
- $\epsilon_{ijk}$ représente l'erreure résiduelle et où : $\forall (i\space ;\space j)\neq (l\space ;\space m) \space : \space \epsilon_{ijk}\perp \epsilon_{lmk}$.              

La décomposition de $\mu_{ij}$ produit un **modèle additif avec interraction**.

En l'occurence, ce modèle est complet.

Remarque : Nombre de paramètres à estimer $1+p+q+p\times q$

### Interaction : 

**Définition :** Variation de l'effet d'un facteur selon l'effet du second.

Soit la moyenne d'un traitement :

$$\mu_{ij}=\mu+\alpha_{i}+\beta_{j}+\gamma_{ij}$$

Si on considère deux niveaux $(i \space ;\space i')$ du même facteur :

$$\forall \space j \space : \space \mu_{ii'}-\mu_{ii'}=\alpha_{ii'}-\alpha_{ii'}+\gamma_{ij}-\gamma_{ii'}$$

On observe que la différence de moyenne entre modélités du même facteur dépend d'un effet du second. Si l'interaction est nulle, le modèle revient à :

$$\mu_{ij}=\mu+\alpha_{i}+\beta_{j}$$

**En l'absence d'interaction, l'écart de moyennes entre modélités du même facteur reste constant à travers les niveaux du second**.

D'un point de vue statistique, cela revient à ce que la **loi commune soit le produit des lois marginales**, avec :

$$P(A) \cap P(B) =P(A)\times P(B)$$

$$\Rightarrow \quad P(A|B)= \frac{P(A) \cap P(B)}{P(B)}=\frac{P(A) \times P(B)}{P(B)}=P(A)$$

$$\Rightarrow \quad P(B|A)= \frac{P(B) \cap P(A)}{P(A)}=\frac{P(B) \times P(A)}{P(A)}=P(B)$$

### L'ANOVA à $2$ facteurs comme modèle de régression : 

Ecriture matricielle :

$$Y=X\theta +\epsilon $$

Avec :

- $X_{(n\times 1+p+q+p\times q)}$ est composée d'une colonne placée à la valeur $1$ (l'intercept), respectivement $p$ et $q$ vecteurs binaires d'appartenance aux facteurs et $p\times q$ vecteurs binaires d'intrices d'interaction.          
- $\theta$ le vecteur ligne des coefficients ($\alpha_{i}$, $\beta_{j}$, $\gamma_{ij}$), de taille $p+q+p\times q$.            
- $\epsilon$ le vecteur ligne de longuer $n$ des erreurs aléatoires.

### Estimation des paramètres du modèle :

Afin de rendre le modèle identifiable, on doit introduire des contraintes d'identification qui peuvent être de deux types :

- Conrtaintes de type **somme** :

$$\sum_{i=0}^{p}\alpha_{i}=0$$

$$\sum_{j=0}^{q}\beta_{i}=0$$

$$\sum_{j=0}^{q}\gamma_{ij}=0$$

$$\sum_{i=0}^{p}\gamma_{ij}=0$$

#### Estimation des paramètres : 

Sous ces contraintes, les paramètres $\{\alpha,\beta,\gamma\}$ sont estimés en minimisant :

$$\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}\left(Y-(\alpha_{i}+\beta_{j}+\gamma_{ij})^{2}\right)$$

Ces estimateurs sont : 

- $\hat{\mu}=Y_{...}$                
- $\hat{\alpha}=Y_{i..}-Y_{...}$               
- $\hat{\beta}=Y_{.j.}-Y_{..}$               
- $\hat{\mu}=Y_{ij.}-Y_{i..}-Y_{.j.}+Y_{...}$                

Remarque : les estimations des paramètres s'interprêtent comme des écarts à la moyenne générale. 

- Contraintes de type **modalité de référence** : en chosissant une modalité de référence pour chaque facteur. 

$$\alpha_{i}=0$$           

$$\beta_{j}=0$$           

$$\forall \space i \space ;\space \gamma_{i1}=0$$           

$$\forall \space j \space ;\space \gamma_{1j}=0$$           

*Dans cet exemple : le traitement $(1\space ;\space 1)$* joue le rôle de référence.

Estimation des paramètres :

Sous ces contraintes, les estimations des paramètres sont données par : 

$$\hat{\mu}=Y_{11.}$$           

$$\hat{\alpha}_{i}=Y_{i1.}-Y_{11.}$$           

$$\hat{\beta}_{i}=Y_{1j.}-Y_{11.}$$           

$$\hat{\gamma}_{ij}=Y_{ij.}-Y_{i1.}-Y_{1j.}+Y_{11.}$$           

Les estimations des paramètres s'interprêtent alors comme des écarts à la moyenne du tratiement de référence. 

Quelles que soient les contraintes choisies, on a : 

$$Y_{ij.}=\hat{\mu}+\hat{\alpha}_{i}+\hat{\beta}_{j}+\hat{\gamma}_{ij}=Y_{ij.}$$

Conséquences :

- Toute variabilité dans les groupes de traitements (intra-groupes) est considérée comme résiduelle.

- Ni l'ajustement, ni la significativité du modèle n'est dépendant des contraintes. 

### Estimation de la variabilité des paramètres :

L'estimateur **sans biais** de la variance est donné par : 

$$S_{n-pq}^{2}=\frac{SCR}{(n-pq)}=\frac{\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{i}}(Y_{ijk}-\hat{Y}_{ijk})^{2}}{(n-pq)}=\frac{\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{i}}(Y_{ijk}-Y_{ij.})^{2}}{(n-pq)}$$

#### Variances des estimateurs :

Soit $r$, le nombre de répétitions dans chaque traitement :

$$\sigma_{\hat{\alpha}_{i}}^{2}=\frac{p-1}{p}\times \frac{\sigma^{2}}{rq}$$

$$\sigma_{\hat{\beta}_{j}}^{2}=\frac{q-1}{q}\times \frac{\sigma^{2}}{rq}$$

$$\sigma_{\hat{\gamma}_{ij}}^{2}=\frac{p-1}{p}\times\frac{q-1}{q}\times \frac{\sigma^{2}}{r}$$

Remarque : dans le cas d'un modèle **équilibré**, $r$ est fixé.

### Décomposition de la variabilité totale :

$$SCT=\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}(Y_{ijk}-Y_{...})^{2}$$

Avec : 

$$SCT=SCR+SCM$$

$$\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}(Y_{ijk}-Y_{...})^{2}=\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}(Y_{ijk}-Y_{ij.})^{2}+\sum_{i=1}^{p}\sum_{j=1}^{q}n_{ij}(Y_{ij.}-Y_{...})^{2}$$

Donc : 

$$SCR=\sum_{i=1}^{p}\sum_{j=1}^{q}\sum_{k=1}^{n_{ij}}(Y_{ijk}-Y_{ij.})^{2}$$

$$SCM=\sum_{i=1}^{p}\sum_{j=1}^{q}n_{ij}(Y_{ij.}-Y_{...})^{2}$$

### Test de significativité globale du modèle :

$\{h_{0} \space : \space Y_ijk=\mu+\epsilon_{ijk}\} \quad \Leftrightarrow \quad \{\forall i \space ;\space \forall j \space : \space \alpha=0 \space ;\space \beta=0 \space ;\space \space ;\space \gamma=0\} \quad \quad$

$\{h_{01} \space : \space Y_ijk=\mu+\alpha_{i}+\beta_{j}+\gamma_{ij}+\epsilon_{ijk}\} \quad \Leftrightarrow \quad \{\exists (i \space ;\space j) \space : \space \alpha \neq 0 \space ;\space \beta \neq 0 \space ;\space \space ;\space \gamma \neq 0\} \quad \quad$

Statistique : 

$$F=\frac{SCM/(pq-1)}{SCR/(n-pq)} \sim \mathcal{F}_{(pq-1 \space ; \space n-pq)}$$

Remarque :

$$SCM \sim_{H_{0}}\chi_{pq-1}^{2}$$

$$SCM \sim_{H_{0}}\chi_{n-pq}^{2}$$

### décomposition de la somme des carrés du modèle :

Hypothèse : plan d'expérience **orthogonal**.

On peut décomposer $SCM$ en :

$$SCM_{A}=\sum_{i=1}^{p}(Y_{i..}-Y_{...})^{2}n_{i.}$$

$$SCM_{B}=\sum_{j=1}^{q}(Y_{.j.}-Y_{...})^{2}n_{.J}$$

$$SCM_{I}=\sum_{i=1}^{p}\sum_{j=1}^{q}(Y_{ij.}-Y_{i..}-Y_{.j.}+Y_{...})^{2}n_{ij}$$

$$\Rightarrow \quad SCM=SCM_{A}+SCM_{B}+SCM_{I}$$

### Test des différents effets :

Hypothèse : Le modèle est globalement significatif.

On teste la significativité des trois composantes de la somme des carrés expliqués :

**Effet du traitement A :**

$H_{0} \space : \{\forall \space i \in [1;...;p] \space : \space \alpha_{i}=0\}$          
$H_{1} \space : \{\exists \space i \in [1;...;p] \space : \space \alpha_{i} \neq 0\}$

Statistique : 

$$F_{A}=\frac{SCM_{A}/(p-1)}{SCR/pq(r-1)}\sim \mathcal{F}_{(p-1\space ; \space pq(r-1))}$$

**Effet du traitement B :**

$H_{0} \space : \{\forall \space i \in [1;...;p] \space : \space \beta_{i}=0\}$          
$H_{1} \space : \{\exists \space i \in [1;...;p] \space : \space \beta_{i} \neq 0\}$

Statistique : 

$$F_{B}=\frac{SCM_{B}/(q-1)}{SCR/pq(r-1)}\sim \mathcal{F}_{(q-1\space ; \space pq(r-1))}$$

**Effet de l'intéraction :**

$H_{0} \space : \{\forall \space (i\ ; \ j) \space : \space \gamma_{ij}=0\}$          
$H_{1} \space : \{\exists \space (i\ ; \ j) \space : \space \gamma_{ij} \neq 0\}$

Statistique : 

$$F_{B}=\frac{SCM_{I}/(p-1)(q-1)}{SCR/pq(r-1)}\sim \mathcal{F}_{(p-1\space ; \space pq(r-1))}$$

### Comparaison des moyennes :

Tests *post-hoc* afin de comparer les moyennes de $Y$ entre des niveaux $i$ et $i'$ d'un facteur pour lesquelles sont différentes :

#### Test de significativité minimale (LSD) de Fisher :

*Dans cet exemple, on s'intéresse au facteur A :*

$H_{0} \space : \space \{Y_{i..}=Y_{i'..}\}$             
$H_{1} \space : \space \{Y_{i..}\neq Y_{i'..}\}$             

Statistique : 

$$\frac{Y_{i..}-Y_{i'..}}{S_{n-p}\sqrt{\frac{1}{n_{i}}+\frac{1}{n_{i'}}}} \sim \mathcal{F}_{(n-p) \space ; \space 1-\alpha/2}$$

Où $S_{n-p}$ est l'estimateur de l'écart type dans un modèle ne comportant que le premier facteur :

$$Y_{ijk}=\mu+\alpha_{i}+\epsilon_{ijk}$$

$$\Rightarrow \quad S_{n-p}=\hat{\sigma}^{2}=\frac{\sum_{i}\sum_{k}(Y_{ijk}-Y_{i..})^{2}}{n-p}$$

### Problème des tests multiples : 

Remarque : plus les niveaux des facteurs sont nomreux et plus important sera le nombre de tests à réaliser.

Pour se prémunir contre l'augmentation du risque $\alpha$ dans le cas de tests multiples, il existe différentes méthodes :

#### Correction de Bonferonni :

Soit $\alpha^{\star}$, la probabilité de faire au moins une erreur de type $1$ sur $l$ tests et $m$ le nombre de tests à réaliser :

On sait que : 

$$\alpha^{\star}\leq\alpha\times l$$
Il faut donc apporter la correction suivante au niveau de significativité attendu : 

$$m=\frac{l(l-1)}{2} \quad \Rightarrow \quad \alpha = \frac{\alpha^{\star}}{m}$$

Remarques : 

- La correction de Bonferroni est très conservative, elle est donc déconseillée en cas de tests successifs nombreux sinon conserver l'hypothèse nulle souvent à tord en raison d'un manque de puissance.

- Des logiciels offrent une *P-Value* ajustée, comme sortie d'un test corrigé par la méthode de Bonferroni ($P^{\star}=P\times l$). 

#### Test de l'étendue de Turkey (HSD) :

Hypothèse : égalité des effectifs dans les groupes.

$$\frac{|Y_{i..}-Y_{i'..}|}{\sqrt{\frac{\hat{\sigma}^{2}}{n/p}}} \sim Q_{(p \space ; \space n-p)} \space ; \space 1-\alpha$$

Où : 

- $\hat{\sigma}^{2}$ est l'estimateur de la variance du modèle.                  
- $n/p$ est l'effectif de chaque groupe (niveaux du facteur).                    
- $Q$ ets la loi de l'étendue studentisée.

Remarque : Plus puissant que la méthode de Bonferroni, le test de Turkey est conseillé lorsque de nombreux tests doivent être réalisés.

### Références
1. Giorgio Russolillo. STA102 : Analyse de la variance à deux facteurs, CNAM.
2. Antoine Godichon-Baggioni. Analyse de la variance à deux facteurs, INSA de Rouen,
http://godichon.perso.math.cnrs.fr/PolyCours5.pdf.
22





```{r}

```

```{r}

```

```{r}

```

```{r}

```

```{r}

```

```{r}

```

