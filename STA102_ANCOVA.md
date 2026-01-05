---
title: "Analyse de la covariance (ANCOVA)"
output: html_notebook
---

Ojectif : expliquer une variable quantitative $Y$ à l'aide de plusieurs autres variables quantitatives et qualitatives, en tenant compte durant l'étude :

- des effets des facteurs.

- des effets de la (les) variable(s) quantitative(s), appelées des **covariables**.

Les comparaisons ne portent plus sur les moyennes au sein de chaque cellule mais sur les paramètres de régression, estimés pour chaque combinaison de niveaux. 

Dans l'exemple suivant on dispose d'un groupe de $n$ individus sur lequels on a mesuré :

- $Y$ : variable à expliquer.                  
- $x$ : variable explicative quantitative.                  
- une variable qualitative appelée le **facteur** et qui présente $p$ niveaux.                  

Remarque : Ajuster le modèle pour $x$ lorsque la covariable est corrélée à $Y$ aide à :

- Réduire la variation résiduelle.                             
- Rendre plus précise la prédiction.                             
- Augmenter la puissance des tests réalisés entre les moyennes du traitement.                    

### Description des données :

- $\forall j\in p \space : \space x_{i} = \{x_{1i}\space ; \space ...\space ; \space x_{n_{i}i}\}$                          
- $\forall i\in p \space : \space y_{i} = \{y_{1i}\space ; \space ...\space ; \space y_{n_{i}i}\}$

$$n=\sum_{i=1}^{p}n_{i}$$

### Modèle sans compensation des effets : 

*Réalisation de $p$ modèles de régression simple de $x$ sur $Y$ en segmentant la population sur chaque niveau du facteur.*

$$Y_{ij}=\mu_{i}+\beta_{i}x_{ij}+\epsilon_{ij}$$

Propriétés : 

- $\forall (i\space ; \space j) \space ; \space \epsilon_{ij} \space iid \space \sim \mathcal{N}(0 \space ; \space \sigma{2})$                         
- $i$ est l'indice du niveau de facteur.                          
- $j$ est l'indice de répétition (numéro de l'observation pour le niveau $i$ du traitement).                          
- Pour $i$ fixé, les estimateurs de $\mu_{i}$ et $\beta_{i}$ sont ceux de la régression linéaire simples : 

$$\forall \space i\in p\space
\begin{cases}
\hat{\beta_{i}}=\frac{\sum_{j=0}^{n_{i}} (x_{i}-\bar{x})(y_{i}-\bar{y})}{\sum_{j=0}^{n_{i}} (x_{i}-\bar{x})^{2}} \\
\hat{\mu}_{i}=Y_{i.}-\hat{\beta}_{i}x_{i.}
\end{cases}
$$

Avec : 

- $Y_{i}$ : la moyenne de la variable à expliquer pour le groupe $i$.
- $x_{i.}$ : la covariable pour le groupe $i$.

En pratique, une première approximation consiste à comparer visuellement les droites de régression pour chaque modèle. Il faut cependant trouver d'autres méthodes pour pousser la comparaison à l'effectif global. 

### Modèle avec décomposition des effets des facteurs :

On décompose entre un **effet global** et d'autres dûs aux différents **niveaux** $i$ du facteur. Soit : 

- $\mu_{i} \space \longrightarrow \space \mu+\alpha_{i}$                    
- $\beta_{i} \space \longrightarrow \space \beta+\gamma_{i}$                    

Avec $\beta$ l'effet global de la covariable et $\gamma$ son effet spécifique au niveau $i$ du facteur. 

Le modèle s'écrit alors :

$$Y_{ij}=\mu+\alpha_{i}+\beta_{i}x_{ij}+\gamma_{i}x_{ij}+\epsilon_{ij} \quad ; \quad \epsilon_{ij} \space iid \sim \quad \mathcal{N}(0 \space ; \space \sigma^{2})$$

Remarque : $\gamma_{ij}$ peut être considéré comme terme d'interaction entre le facteur et la covariable. 

#### Ecriture matricielle :

$$y=X\theta+\epsilon$$
Soit :

```math
\mathbf y
=
\begin{pmatrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{pmatrix}
\qquad
\mathbf X
=
\begin{pmatrix}
1 & 0 & \cdots & 1 & x_{1} & 0 & \cdots & Y_{1p} \\
1 & 1 & \cdots & 0 & x_{2} & Y_{2,1} & \cdots & 0  \\
\vdots & \vdots &  & \vdots & \vdots & \vdots &  & \vdots \\
1 & 0 & \cdots & 1 & x_{n_i} & 0 & \cdots &  Y_{n_i p} \\
\end{pmatrix}
\qquad
\boldsymbol\theta
=
\begin{pmatrix}
\mu\\
\alpha_1\\
\vdots\\
\alpha_p\\
\beta\\
\gamma_1\\
\vdots\\
\gamma_p
\end{pmatrix}
\qquad
\boldsymbol\varepsilon
=
\begin{pmatrix}
\varepsilon_1\\
\varepsilon_2\\
\vdots\\
\varepsilon_n
\end{pmatrix}
```


Avec $X$, une X matrice $n × (1 + p + q + p \times q)$ et $theta$ le vecteur ligne des coefficients ($\mu,\alpha_{i},\beta,\gamma_{i}$).

#### Contraintes d'identifibilité du modèle :

L'estimation des paramètres ne peut être réalisée qu'avec l'ajout de **deux contraintes** :

- Contrainte de type **somme** :

$$\sum_{i=1}^{p}n_{i}\alpha_{i}=\sum_{j=1}^{q}n_{i}\beta=0$$

- Contrainte de type **niveau de référence** :

$$\alpha_{1}=\beta_{1}=0$$

Remarque : Sous R, le niveau par défaut est le premier.

Les modèles avec et sans décomposition expriment de façon différente les paramètres. Soit, en prenant l'exemple d'une covariable $1$  :

$\forall \space i \in [1 \space ;...;\space p]$ :

- $\mu=\mu_{1}$                                 
- $\beta=\beta_{1}$                                 
- $\alpha=(\mu_{i}-\mu_{1})$                                 
- $\gamma=(\beta_{i}-\beta_{1})$                                 

#### Prévisions :

*Les prévisions sont indépendantes des contraintes et correspondent à celles de $p$ régressions linéaires.*

#### Estimateur de la variance des résidus : 

$$S_{n-2p}^{2}=\frac{SCR}{n-2p}=\frac{\sum_{i=1}^{p}\sum_{j=1}^{n_{i}}(Y_{ij}-\hat{Y}_{ij})^{2}}{n-2p}$$

On note sa réalisation $\hat{\sigma}^{2}$.

### Décomposition de la variabilité totale :

$$SCM=SCM+SCR$$






















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

```{r}

```

