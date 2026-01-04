---
title: "ANOVA à un facteur"
output: html_notebook
---

### Contexte :

Etude de l'influence d'un ou plusieurs facteurs sur la moyenne $\mu$ d'une variable quantitative. Deux optiques possibles : 

- Expliquer une variable quantitative $Y$ au moyen d'un facteur à $K$ modélités, appelés des niveaux.

- Comparer des populations ou différentes conditions expérimentales.

### Données : 

Soit $K$ échantillons de taille $\{n_{1},...,n_{k}\}$ (cf. partition de la population totale en $k$ groupes) avec :

- $\sum_{i=1}^{k}n_{i}=n$

- La moyenne de $Y$ pour le groupe $i$ : $Y_{i.}=\frac{1}{n_{i}}\sum_{j=1}^{n_{i}}y_{ij}$

- La moyenne de $Y$ pour l'ensemble des observations : $Y_{..}=\frac{1}{n}\sum_{i=1}^{K}\sum_{j=1}^{n_{i}}y_{ij}$

### Modélisation : 

$$Y_{ij}=\mu_{i} + \epsilon_{ij} \quad ; \quad i\in[1,...,k], j\in[1,...,n_{i}]$$

Avec : $\mu_{i}$ la **moyenne théorique** (ou **effet**) de la variable à expliquer dans le groupe $i$.

Seconde modélisation :

On décompose $\mu$ afin de mettre en évidence un éventuel effet global des effets marginaux :

$$\mu_{i}=\mu + \alpha_{i} \quad ; \quad i\in[1,...,k]$$

Le modèle devient alors : 

$$Y_{ij}=\mu + \alpha_{i} + \epsilon_{ij} \quad ; \quad i\in[1,...,k], j\in[1,...,n_{i}]$$

### Hypothèses sur les erreurs :

Dans les deux cas, on fait les hypothèses suivantes :

- $\forall i \ : \epsilon \quad  iid$          
- $\forall i,\ \epsilon_{i}\sim \mathcal{N}(0,\sigma^{2})$

Donc : 

- $\forall \{i\ ; \ j\}$ : les $y_{ij}$ sont des variables aléatoires centrées sur $\mu_{i}$.       
- $\forall \{i\ ; \ j\} \ : \ y_{ij} \quad iid$

### Estimation des paramètres :

#### Premier modèle :

$K$ paramètres à estimer :

$$
\hat{\mu}_{i}=\frac{1}{n_{i}}\sum^{n_{i}}_{j=1}y_{ij}
$$

#### Validation de l'estimateur par la méthode des moindres carrés :

##### Ecriture matricielle : 

$$
Y=\mu_{1}x_{1}+...+\mu_{k}x_{k}+\epsilon
$$

où $x_{i}$ est une indicatrice d'appartenance au groupe $i$, telle que :

$$
x_{ij}=\begin{cases}
1 \quad \rightarrow \quad x_{ij} \in i \\
0 \quad \rightarrow \quad x_{ij} \not\in i
\end{cases}
$$

Remarque : ce modèle **peut être considéré comme un model de régression linéaire multiple** :

- $k$ variables explicatives         
- Les $\mu_{i}$ correspondent aux $\beta_{i}$         
- L'**intercept** est nulle         
- $\epsilon$ est le terme d'erreur         

Forme  matricielle : 

$$
Y =
\begin{pmatrix}
Y_{11} \\ 
\vdots \\ 
Y_{1 n_1} \\ 
Y_{21} \\ 
\vdots \\ 
Y_{2 n_2} \\ 
\vdots \\ 
Y_{k1} \\ 
\vdots \\ 
Y_{k n_k}
\end{pmatrix}
\quad
X =
\begin{pmatrix}
1 & 0 & 0 & \cdots & 0 \\
\vdots & \vdots & \vdots &        & \vdots \\
1 & 0 & 0 & \cdots & 0 \\
0 & 1 & 0 & \cdots & 0 \\
\vdots & \vdots & \vdots &        & \vdots \\
0 & 1 & 0 & \cdots & 0 \\
\vdots &        &        & \ddots & \vdots \\
0 & 0 & 0 & \cdots & 1 \\
\vdots & \vdots & \vdots &        & \vdots \\
0 & 0 & 0 & \cdots & 1
\end{pmatrix}
\quad
\beta =
\begin{pmatrix}
\mu_1 \\ 
\mu_2 \\ 
\vdots \\ 
\mu_k
\end{pmatrix}
\quad
\varepsilon =
\begin{pmatrix}
\varepsilon_{11} \\ 
\vdots \\ 
\varepsilon_{1 n_1} \\ 
\varepsilon_{21} \\ 
\vdots \\ 
\varepsilon_{2 n_2} \\ 
\vdots \\ 
\varepsilon_{k1} \\ 
\vdots \\ 
\varepsilon_{k n_k}
\end{pmatrix}
$$

#### Estimation par la méthode des moindres carrés : 

$$b=X(X^{T}X)^{-1}X^{T}Y$$

Avec : 

$$
X^T X =
\begin{pmatrix}
n_1 & 0 & \cdots & 0 \\
0 & n_2 & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & n_k
\end{pmatrix}
\quad
X^T Y =
\begin{pmatrix}
n_1 \bar{Y}_1 \\
n_2 \bar{Y}_2 \\
\vdots \\
n_k \bar{Y}_k
\end{pmatrix}
\quad
(X^T X)^{-1} =
\begin{pmatrix}
\frac{1}{n_1} & 0 & \cdots & 0 \\
0 & \frac{1}{n_2} & \cdots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \cdots & \frac{1}{n_k}
\end{pmatrix}
$$

Donc : 

$$
\hat{\mu}_{i}=Y_{i.}=\frac{1}{n_{i}}\sum^{n_{i}}_{j=1}y_{ij}
$$

#### Second modèle : 

$(k+1)$ paramètres à estimer ($\mu,\alpha_{1},...,\alpha_{k}$) :

Remarque : il existe une infinité de décompositions possibles. Il faut donc choisir des **conditions d'identifiabilité** sur les paramètres.

##### Ecriture matricielle : 

$$Y_{ij} = \mu+\alpha_{1}x_{1}+...+\alpha_{1}x_{1}+\epsilon_{ij}$$

Avec : 

- $\mu$ l'ordonnée à l'origine.
- $k$ variables explicatives $\{x_{1},...,x_{k}\}$.
- Les $\alpha_{i}$ sont identifiés aux coefficients de régression. 

Forme matricielle :

$$
Y=X\beta+\epsilon
$$

$$
Y =
\begin{pmatrix}
Y_{11} \\ 
\vdots \\ 
Y_{1 n_1} \\ 
Y_{21} \\ 
\vdots \\ 
Y_{2 n_2} \\ 
\vdots \\ 
Y_{k1} \\ 
\vdots \\ 
Y_{k n_k}
\end{pmatrix}
\quad
X =
\begin{pmatrix}
1 & 1 & 0 & 0 & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots &        & \vdots \\
1 & 1 & 0 & 0 & \cdots & 0 \\
1 & 0 & 1 & 0 & \cdots & 0 \\
\vdots & \vdots & \vdots & \vdots &        & \vdots \\
0 & 0 & 0 & 0 & \cdots & 1 \\
\vdots & \vdots & \vdots & \vdots &        & \vdots \\
1 & 0 & 0 & 0 & \cdots & 1 \\
\vdots & \vdots & \vdots & \vdots &        & \vdots \\
1 & 0 & 0 & 0 & \cdots & 1
\end{pmatrix}
\quad
\beta =
\begin{pmatrix}
\mu \\ 
\alpha_1 \\ 
\alpha_2 \\ 
\vdots \\ 
\alpha_k
\end{pmatrix}
\quad
\varepsilon =
\begin{pmatrix}
\varepsilon_{11} \\ 
\vdots \\ 
\varepsilon_{1 n_1} \\ 
\varepsilon_{21} \\ 
\vdots \\ 
\varepsilon_{2 n_2} \\ 
\vdots \\ 
\varepsilon_{k1} \\ 
\vdots \\ 
\varepsilon_{k n_k}
\end{pmatrix}
$$

#### Estimation par la méthode des moindres carrés : 

Remarque : dans ce modèle on a un problème d’identification, car $rang(X) = k < k +1$ (la première colonne de X étant égale à la somme des autres colonnes) lorsqu’on à (k + 1) paramètres inconnus à estimer.

Les **contraintes** que nous devons appliquer sont de deux types : 

- Contraintes appliquées à la somme des effets : 

$$
\sum^{k}_{i=1}\alpha_{i}=0 \quad \Rightarrow \quad \hat{\mu}=\frac{1}{k}\sum^{k}_{i=1}Y_{i.}
$$
$$
\hat{\alpha}_{i}=Y_{i.}-\hat{\mu}
$$

**Deux cas possibles :**

- **dispositif équilibré :** $\forall \{i\ ; \ j\}\in k \ ; \ i \neq j \ : \ n_{i}= n_{j}$ :

$$
\hat{\mu}=\frac{1}{k}\sum^{k}_{i=1}Y_{i.}=Y_{..}
$$

$$
\hat{\alpha}_{i}=Y_{i.}-\hat{\mu}=Y_{i.}-Y_{..}
$$

- **dispositif déséquilibré :** $\forall \{i\ ; \ j\}\in k \ : \ i \neq j \Rightarrow n_{i}\neq n_{j}$ :

La contrainte devient : $\sum^{k}_{j=1}n_{i}\alpha_{i}=0$

- Contrainte appliquée à une modalité de référence : on estime les effets par rapport à une modalité de référence (notée ici $k$) : 

$$
\alpha_{k}=0
$$

$$
\hat{\mu}_{k}=Y_{k.}
$$

$$
\hat{\alpha}_{i}=Y_{i.}-Y_{k.}
$$

### Prédiction : 

L'estimation de l'espérance mathématique de $Y$ (moyenne selon le groupe d'appartenance) **ne dépend pas du modèle**. Dans tous les cas : 

$$
\hat{Y}_{ij}=Y_{i.}
$$

#### Démonstration : 

Modèle 1 :

$$
Y_{ij}=\mu_{i}+\epsilon_{ij}
$$

$$
\hat{Y}_{ij}=\hat{\mu}_{i}=Y_{i.}
$$

Modèle 2 :

- Contrainte : $\sum^{k}_{i=1}n_{i}\alpha_{i}=0$

$$
Y=\mu+\alpha_{i}+\epsilon_{ij}
$$

$$
\hat{Y}_{ij}=\hat{\mu}+\hat{\alpha}_{i}=\hat{\mu}+Y_{i.}-\hat{\mu}=Y_{i.}
$$

- Contrainte : $\alpha_{k}=0$

$$
Y_{ij}=\mu+\alpha_{i}+\epsilon_{ij}
$$

$$
\hat{Y}_{ij}=\hat{\mu}+\hat{\alpha}_{i}=Y_{k\ .}+Y_{i.}-Y_{k.}=Y_{i.}
$$

### Estimateur de \sigma^{2} :

Rappel : selon le modèle d'*ANOVA*, toutes les observations d'un groupe devrait approcher la même valeur. La variabilité intra-groupe est considérée comme une erreur aléatoire :

Conséquences : 

- Résidu aléatoire pour $Y_{ij}$ : $Y_{ij}-\hat{Y}_{ij}=Y_{ij}-Y_{i.}$

- Estimateur sans biais de la variance résiduelle (\sigma^{2}) : 

$$
S^{2}_{n-k}=\frac{1}{n-k}\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Y_{ij}-Y_{i.})^{2}
$$

### Décomposition de la variance totale : 

Remarque : le modèle d'*ANOVA* étant un modèle de régression, on peut décomposer la variance entre une partie expliquée par le modèle et une autre résiduelle. 

$$
SCT=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Y_{ij}-Y_{..})^{2}
$$

$$
SCR=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Y_{ij}-Y_{i.})^{2}
$$

$$
SCF=\sum^{k}_{i=1}n_{i}(Y_{i.}-Y_{...})^{2}
$$

On retrouve : $SCT=SCR+SCF$

$$
SCT=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{..}\right)^{2}=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}+Y_{i.}-Y_{..}\right)^{2}
$$

$$
=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)^{2}+\left(Y_{i.}-Y_{..}\right)^{2}+2\left(Y_{ij}-Y_{i.}\right)\left(Y_{i.}-Y_{..}\right)
$$

$$
=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)^{2}+\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{i.}-Y_{..}\right)^{2}+2\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)\left(Y_{i.}-Y_{..}\right)
$$

$$
=\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)^{2}+\sum^{k}_{i=1}n_{i}\left(Y_{i.}-Y_{..}\right)^{2}+2\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)\left(Y_{i.}-Y_{..}\right)
$$

$$
=SCR+SCF+2\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)\left(Y_{i.}-Y_{..}\right)
$$

Avec :

$$
2\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}\left(Y_{ij}-Y_{i.}\right)\left(Y_{i.}-Y_{..}\right)=0
$$

$$
\Rightarrow \quad SCT=SCR+SCF
$$

### Test de significativité globale du modèle : 

Test de l'effet du facteur sur $Y$.

$H_{0} \ : \ \forall i\neq j\in k \ ; \ y_{i} = y_{j} \quad \Leftrightarrow \quad \{Y_{ij}=\mu+\alpha_{i}\} \quad \Leftrightarrow \quad \{\forall i \ ; \ \alpha_{i}=0\}$
$H_{1} \ : \ \exists i,i'\in k \ ; \mu_{i} \neq \mu_{j} \quad \Leftrightarrow \quad \{Y_{ij}=\mu+\alpha_{i}+\epsilon_{ij}\}\quad \Leftrightarrow \quad \{\exists \alpha_{i}\neq0\}$ 

Statistique : 

$$
F=\frac{SCE/(k-1)}{SCR/(n-1)}\sim \mathcal{F}_{(k-1\ ; \ n-1)}
$$

Car, sous $H_{0}$ :

$$
\frac{\sum^{n_{i}}_{j=1}(Y_{ij}-Y_{i.})^{2}}{\sigma^{2}} \sim \chi^{2}_{n_{i}-1}
$$

$$
\Longrightarrow \quad\frac{SCR}{\sigma{2}}=\frac{\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Y_{ij}-Y_{i.})^{2}}{\sigma^{2}} \sim \chi^{2}_{n-k}
$$

$$
\frac{SCT}{\sigma{2}}=\frac{\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Y_{ij}-Y_{..})^{2}}{\sigma^{2}} \sim \chi^{2}_{n-1}
$$

$$
\frac{SCF}{\sigma{2}}=\frac{\sum^{k}_{i=1}n_{i}(Y_{i.}-Y_{..})^{2}}{\sigma^{2}} \sim \chi^{2}_{n-1}
$$

Remarque : rejetter $H_{0}$ signifie admettre que $X$ joue une rôle significatif sur $Y$.

#### Tableau d'analyse de la variance : 

| Source de variation | Degrés de liberté (DF) | Somme des carrés | Carré moyen | Statistique F |
|--------------------|------------------------|-----------------|-------------|---------------|
| Modèle             | $k-1$                  | $\mathrm{SCF}$  | $\dfrac{\mathrm{SCF}}{k-1}$ | $F = \dfrac{\mathrm{SCF}/(k-1)}{\mathrm{SCR}/(n-k)}$ |
| Erreur             | $n-k$                  | $\mathrm{SCR}$  | $\dfrac{\mathrm{SCR}}{n-k}$ |               |
| Totale             | $n-1$                  | $\mathrm{SCT}$  |             |               |

### Validation des hypothèses du modèle :

Remarque : l'*ANOVA* est réputée assez robuste aux écarts à ses hypothèses.

- Normalité : principalement visuelle ou par tests d'hypothèse (à documenter). L'*ANOVA* robiste à une violation modérée de normalité.

- Homoscédasticité : le test classique est celui de Barlett (à documenter). L'*ANOVA* est par contre **assez sensible à une violation de cette hyptohèse** ; il est donc nécessaire de la valider systématiquement. 

**Attention :** le test de Barlett est **sensible à la non-normalité**. En cas de doute, lui préférer celui de Levene :

Où les $Z_{i}$ sont des mesures de la distance des observations à leurs moyennes : 

$$
W=\frac{\frac{\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Z_{i.}-Z_{..})^{2}}{k-1}}{\frac{\sum^{k}_{i=1}\sum^{n_{i}}_{j=1}(Z_{ij}-Z_{i.})^{2}}{n-k}} \sim \mathcal{F}_{(k-1 \ ; \ n-k)}
$$

Ce test suppose l'**homoscédasticité**. Il existe plusieurs variantes : 

- distances aux moyennes : $Z_{ij}=(Y_{ij}-Y_{i.})^{2}$

- distances en valeur absolue : $Z_{ij}=|Y_{ij}-Y_{i.}|$

- distances à la médiane : $Z_{ij}=(Y_{ij}-\tilde{Y}_{i.})^{2}$

### ANOVA non paramétrique : 

A documenter 

### Références

1. Giorgio Russolillo. STA102 : Analyse de la variance à un facteur, le CNAM.

2. Marie Chavent. Analyse de la variance, Université de Bordeaux, http://www.math.
u-bordeaux.fr/~mchave100p/wordpress/wp-content/uploads/2013/10/CoursModStat_
C6.pdf.

3. E. Lebarbier, S. Robin (2004). Exemples d’application du mod`ele linéaire, Institut
National Agronomique Paris - Grignon http://www2.agroparistech.fr/IMG/pdf/
ModLin-exemples.pdf

```{r}

```

```{r}

```

```{r}

```

```{r}

```
