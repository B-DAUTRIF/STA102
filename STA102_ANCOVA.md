---
title: "Analyse de la covariance (ANCOVA)"
output: html_notebook
---

Ojectif : expliquer une variable quantitative $Y$ à l'aide de plusieurs autres variables quantitatives et qualitatives, en tenant compte :

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

$$
\mathbf{X}
=
\begin{pmatrix}
1 & 0 & \cdots & 1 & x_1 & 0 & \cdots & Y_{1p} \\
1 & 1 & \cdots & 0 & x_2 & Y_{21} & \cdots & 0 \\
\vdots & \vdots &  & \vdots & \vdots & \vdots &  & \vdots \\
1 & 0 & \cdots & 1 & x_n & 0 & \cdots & Y_{np}
\end{pmatrix}
$$



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

$$SCT=\sum_{i=1}^{p}\sum_{j=1}^{n_{i}}(Y_{ij}-Y_{..})^{2}$$

$$SCR=\sum_{i=1}^{p}\sum_{j=1}^{n_{i}}(Y_{ij}-\hat{Y}_{ij})^{2}$$

$$SCM=\sum_{i=1}^{p}\sum_{j=1}^{n_{i}}(\hat{Y}_{ij}-Y_{..})^{2}$$

On produit le tableau d'analyse de la variance :

| Source  | Degrés de liberté | Somme des carrés | Carrés moyens |
|--------|------------------|-----------------|---------------|
| Modèle | $2p - 1$ | $\displaystyle \mathrm{SCM} = \sum_{i=1}^{p}\sum_{j=1}^{n_i}\left(\hat Y_{ij} - Y_{..}\right)^2$ | $\displaystyle \frac{\mathrm{SCM}}{2p-1}$ |
| Résidu | $n - 2p$ | $\displaystyle \mathrm{SCR} = \sum_{i=1}^{p}\sum_{j=1}^{n_i}\left(Y_{ij} - \hat Y_{ij}\right)^2$ | $\displaystyle \frac{\mathrm{SCR}}{n-2p}$ |
| Total  | $n - 1$ | $\displaystyle \mathrm{SCT} = \sum_{i=1}^{p}\sum_{j=1}^{n_i}\left(Y_{ij} - Y_{..}\right)^2$ | $\displaystyle \frac{\mathrm{SCT}}{n-1}$ |

### Test de significativité du modèle complet : 

- $H_{0} \space : \space \{Y_{ij}=\mu+\epsilon_{ij}\}$

- $H_{1} \space : \space \{Y_{ij}=\mu+\alpha_{i}+\beta x_{ij}+\gamma_{i} x_{ij}+\epsilon_{ij}\}$

$$F=\frac{SCM/n-(2p-1)}{SCR/(n-2p)} \sim \mathcal{F}_{(2p-1 \space ; \space n-2p)}$$

### Sous-modèles du modèle complet : 

Les modèles ci-dessous sont **successivement imbriqués** :

#### Modèle nul :

Modèle constant où aucun effet n'est présent :

$$Y_{ij}=\mu+\epsilon_{ij}$$

#### Modèle A :

Uniquement l'effet du facteur y est présent :

$$Y_{ij}=\mu+\alpha_{i}+\epsilon_{ij}$$

Remarque : Le modèle A est équivalent à une ANOVA à 1 facteur 

#### Modèle B :

Aucun effet du facteur n'est présent. 

$$Y_{ij}=\mu+\beta x_{ij}+\epsilon_{ij}$$

Remarque : Ce modèle est équivalent à une régression simple. 

#### Modèle AB :

On observe un effet à la fois de la covariable et du traitement : les moyennes de $Y$ dans les différentes modalités du traitement sont différentes. Cette différence reste **constante** à travers les valeurs de la covariable : **c'est le modèle de l'ANCOVA**.

$$Y_{ij}=\mu+\alpha_{i}+\beta x_{ij}+\epsilon_{ij}$$

### Modèle C :

Où on observe une différence de l'effet du traitement selon les modalités de la covariable (et vice-versa) :

$$Y_{ij}=\mu+\alpha_{i}+\beta x_{ij}+\gamma_{i} x_{ij}+\epsilon_{ij}$$

### Décomposition de la somme des carrés du modèle :

Hypothèse : significativité globale du modèle (invalidation du modèle nul).

**Attention :** Dans un plan d'expérience non orthogonal, il est impossible en l'état de décomposer la variabilité du modèle.

Plusieurs méthodes peuvent être utilisées pour obtenir les sommes des carrés :

**Définition :** La **réduction** est la mesure de la différence observée dans le cas de modèles imbriqués lors du passage de l'un modèle à l'autre. 

$$R(M_{1}|M_{2})=SCR_{M_{1}}-SCR_{M_{2}}$$

#### Somme des carrés de type I :

Construction du modèle complet par ajout sucessifs d'effets. La séquence des effets ajoutés peut être adaptée mais conditionne le modèle final, mais leur somme est égale à $SCM$.

#### Somme des carrés de type III :

Comparaison du modèle complet à d'autres construits par éviction successive d'un paramètre. La somme des carrés de type III est indépendante de l'ordre mais n'est pas égale à $SCM$.

### Tests des effets :

Soit le modèle de référence $M_{ref}$ et deux sous-modèles $M_{1}$ et $M_{2}$ tel que :

$$M_{1} \subset M_{2} \subset M_{ref}$$

- $H_{0} \space : \space \{SCR_{M_{1}}-SCR_{M_{2}}=0\}$

- $H_{1} \space : \space \{SCR_{M_{1}}-SCR_{M_{2}}\neq 0\}$

Statistique :

$$\frac{(SCR_{M_{1}}-SCR_{M_{2}})/(df_{M_{1}}-df_{M_{2}})}{SCR_{M_{ref}}/df_{M_{ref}}} \sim \mathcal{F}_{(df_{M_{1}}-df_{M_{2}}\space ;\space df_{M_{ref}})}$$

Remarque : Le dénominateur correspond à l'estimation de la variance résiduelle du modèle de référence. 

#### Aspects pratiques : 

On évalue généralement d'abord l'interaction (avec la coïncidence des somme de carrés de type I et III) : 

- si l'interaction est significative on a trouvé le meilleur modèle.

- sinon on teste les effets du facteur et de la coviriable indépendamment avant de les conserver.

### Comparaison des moyennes :

Evaluation de la significativité statistique dans les $p$ niveaux du facteur. 

- Moyenne **non ajustée** : dépend de la valeur de la covariable.

$$\hat{\mu}_{i.}=\hat{\mu}+\hat{\alpha}_{i}+\hat{\beta}x_{ij}+\hat{\gamma}_{i}x_{i.}$$

Remarque : cette valeur est la même que celle d'un groupe présentant une valeur moyenne de covariable dans le $i^{eme}$ groupe.

- Moyenne **ajusté** : réalise une comparaison des effets du facteur **à parité de la covariable**.

$$\tilde{Y}_{i.}=\hat{\mu}+\hat{\alpha}_{i}+\hat{\beta}x_{ij}+\hat{\gamma}_{i}x_{..}$$

Remarque : Cette approche s'avère risquée si l'interaction est significative. Il faut dans ce cas utiliser plusieurs valeurs de référence.

### Références :

1. Giorgio Russolillo. STA102 : Analyse de la variance à deux facteurs, CNAM.

2. E. Lebarbier, S. Robin (2004). Exemples d’application du modèle linéaire Institut
National Agronomique Paris – Grignon.



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

