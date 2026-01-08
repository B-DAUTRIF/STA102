---
title: 'Regression_multiple'
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

library(prettyR)
library(car)
library(ggplot2)
library(ggfortify)
library(lmtest)
library(olsrr)
library(corrplot)
```

## Données : 

```{r Import}
# Import du jeu de donnees :
df <- read.csv2(
  "D:/Etudes/2025_2026/STA102/3.Regression multiple/TP_regression_multiple/ozone.csv"
)
head(df)
```


```{r Synthese du jeu de donnes}
# Affichage d'informations statistiques generales sur le jeu de donnees :
describe(df)
```

### Corrélations : 

```{r Recherche de correlations dans les donnes}
# Formule du calcul du coefficient de correlation r :
# sqrt(cov(y,x)/(var(x)*var(x)))

# Representation graphique des correlations entre variables :
pairs(df)
# Matrice de correlations :
cor_df=cor(df)
cor_df
# Matrice graphique de correlations :
corrplot(cor_df, order = 'hclust', method = 'ellipse', type = 'lower')
```

## Modèle de régression multiple : 

$$\mathcal{y}_{i}=\beta_{0}+\sum_{j=1}^{p}\beta_{j}x_{i}+e_{i}$$

Soit (en écriture matricielle) :

$$Y=X\beta+\epsilon$$

Avec :

$$
y=
\begin{pmatrix}
y_1 \\
y_2 \\
\vdots \\
y_n
\end{pmatrix}_{n\times 1}
\quad
X=
\begin{pmatrix}
1 & x_{11} & \cdots & x_{1p} \\
1 & x_{21} & \cdots & x_{2p} \\
\vdots & \vdots &        & \vdots \\
1 & x_{n1} & \cdots & x_{np}
\end{pmatrix}_{n\times p}
\quad
\beta=
\begin{pmatrix}
\beta_0 \\
\beta_1 \\
\vdots \\
\beta_p
\end{pmatrix}_{p\times 1}
\quad
\varepsilon=
\begin{pmatrix}
\varepsilon_1 \\
\varepsilon_2 \\
\vdots \\
\varepsilon_n
\end{pmatrix}_{n\times 1}
$$

```{r Creation du modele lineaire}
# Creation du modèle complet incorporant toutes les variables :
modele=lm(maxO3~.,data=df)
```

### Création des matrices caractéristiques :

```{r Matrices caracteristiques}
# Creation de la matrice X des observations :
X=cbind(1,as.matrix(df[,-1]))
# Creation de la matrice Y des resultats :
Y=as.matrix(df[,1])
# Creation de la matrice A de l'EMC :
A=solve(t(X)%*%X)%*%t(X)
# Creation de hat-matrice H :
H=X%*%solve(t(X)%*%X)%*%t(X)
# Creation de la matrice M :
M=diag(dim(H)[1])-H
```

### Synthèse du modèle : 

#### Paramètres du modèle : 

```{r Fonctions de synthese du modele}
# Vecteur de parametre beta par la methode du calcul matriciel :
beta=A%*%Y
# Verification de la validite du calcul :
# le vecteur obtenu de cette facon est une approximation a 10 decimales du vecteur 
# fourni par le modele 
if(all(modele$coefficients-beta<10**-10)){print("Calcul correct")}
```

```{r Fonctions de synthese du modele}
# Synthese :
summary(modele)
# Coefficients avec IC :
Confint(modele)
# Table des sommes de carres :
anova(modele)
```
Remarque : 

- Si la **P-value** affiché par la commande `summary` est inférieure à la limite visée pour le risque de première espèce ($\alpha$), on en déduit qu'au moins une variable explicative a un effet significatif sur $Y$.

- Dans le modèle complet, l'absence de significativité d'un coefficient peut signifier soit la non significativité de l'effet de la variable, soit une redondance d'information.

### Corrélation des vecteurs colonnes de $X$ avec le vecteur réponse $Y$ :

```{r Recherche de correlations dans le modele}
# Creation de la matrice de correlation entre variables explicaives et Y :
cor_beta_Y=cor(df[,-1],df[,1])
corrplot(cor_beta_Y)
```
- Une correlation forte ($<0.5$) oriente le défaut de significativité du modèle vers la structure des corrélations entre variables explicatives. Il faut alors faire retour au *corrélogramme* en début d'analyse.

## Recherche de multicolinéarité :

On fait le diagnostique afin d'avoir la mesure de la multicolinéarité pour chaque variable explicative. On utilise pour cela les quantités suivantes : 

### Tolérance : 

$$Tol(x_{p})=1-R_{p}^{2}$$

où $R_{p}^{2}$ est le coefficient d'une variable explicative sur les autres.

Propriétés :

- $Tol(x_{p})\in[0\space ; \space 1]$

- Pour $Tol(x_{p})=0$, la variable est **combinaison linéaire des autres**. 

### Facteur d'Inflation de la Variance :

$$VIF=\frac{1}{Tol(x_{p})}$$

Propriétés : 

- $VIF \in [1\space ;+\infty]$
- $Tol=1 \ \Rightarrow \ VIF=1$

**Rappel :** La variance des paramètres peut être déduite par : 

$$\hat\sigma^{2}_{B_{p}}=\frac{\hat{\sigma^{2}}}{\Sigma^{n}_{i=1}(x_{ip}-\bar{x}_{p})^{2}} \times VIF(x_{p})$$

Ordres de grandeur : il est généralement souhaitable de respecter :

- $Tol \geq 0.5$
- $VIF_{x_{p}} \leq 2$

```{r Recherche de multicolinearite}
# Affichage de la Tolerance et du Facteur d'Inflation de la Variance (VIF) :
vif_tol_modele = ols_vif_tol(modele)
vif_tol_modele
# Autres mesures de colinearite :
round(ols_eigen_cindex(modele), 2)

# Selection des variables suspectes : 
# - En terme de Tolerance :
Liste_tol=vif_tol_modele[which(vif_tol_modele$Tolerance<=0.5),"Variables"]
# - En terme de VIF :
Liste_VIF=vif_tol_modele[which(vif_tol_modele$VIF>=2),"Variables"]
# ATTENTION : les deux ecritures sont equivalentes et a utiliser suivant le contexte.
```

### ACP des vecteurs de variables explicatives :

On peut analyser le sous-espace vectoriel des variables explicatives. Pour cela on décrit $(X^{T}X)^{-1}$ en terme de ses **composantes principales** : 

- *Condition Index* (*Condition Number*) : racine du rapport entre la plus grande et la plus petite des composantes (valeur propre)  dans $(X^{T}X)^{-1}$

Remarque : $CI>30$ suppose un risque de multicolinéarité.

```{r ACP espace de variables explicatives}
# Caclul des valeurs propres des vecteurs colonnes de X (variables explicatives) :
eigen_cindex_modele = round(ols_eigen_cindex(modele), 2)
# Selection des variables problematiques dans la liste du modele: 
Liste_CI = labels(modele$terms)[index(eigen_cindex_modele[which(eigen_cindex_modele$`Condition Index`>=30),])]
```

## Sélection des variables :

### Best Subset Regression : méthode automatisée **ascendante exhaustive** de sélection des variables, par comparaison systématique des différents modèles emboités parcourus par l'algrithme. 

```{r Methode ascendante exhaustive de selection des variables}
# Sélection exhaustive ascendante des variables :
ols_step_best_subset(modele)
# Graphiques du processus de selection automatique :
bsr = ols_step_best_subset(modele)
plot(bsr)
```

### Best Subset Regression : méthode **pas-à-pas  descendante exhaustive**, avec ici aussi une comparaison systématique des différents modèles parcourus par l'algrithme. 

```{r Methode pas a pas descendante de selection des variables}
# Sélection pas a pas descendante de variables :
sba=ols_step_backward_aic(model = modele, detail = T)
plot(sba)
# On compare le resultat au modele complet en terme de tolerance et de VIF :
modele_r=sba$model
ols_vif_tol(modele_r)
```



```{r }
# Représentation graphique du processus de sélection : 
bsr = ols_step_best_subset(modele)
plot(bsr)

# Pour l'exercice, une autre methode de selection : descendante, pas a pas,
# basee sur l'AIC :
# temp = ols_step_backward_aic(modele, details = T,progress = T)
```

```{r Evaluation du modele restreint}
# Synthese du modele restreint :
summary(modele_r)
# Sommes des carres et analyse de variance dans le modele restreint :
anova(modele_r)
# Etude des IC : 
confint(modele_r)
```

On remarque qu'effectivement aucun des IC des coefficients ne contient 0.

## Diagnostique des résultats :

### Analyse des résidus : 

Affichage des graphiques de base pour les résidus :
```{r Graphiques standards des residus du modele}
# Affichage des graphiques :
# Residus VS predictions :
autoplot(modele_r, which = 1, ncol = 2, label.size = 3, label.hjust = -0.8, label.n = 4)+
# Q-Q plot :
autoplot(modele_r, which = 2, ncol = 2, label.size = 3, label.hjust = -0.8, label.n = 4)

# Lot de graphiques standards des residus du modele :
plot(modele_r)
# ATTENTION : commande a n'utiliser qu'en console.
```

### Affichage de la synthese des residus par observations :

```{r Synthese des residus}
# Synthese detaille par observation des residus : 
broom::augment(modele_r)
```

### Tests statistiques des residus :

On exécute des tests statistiques sur les résidus pour vérifier les hypothèses du modèle :

- Recherche d'autocorrélation avec le test de Durbin-Watson ($H_0:$ Pas d'autocorrélation).

- Test d'homoscédasticité de Breusch-Pagan ($H_0:$ Homoscédasticité).

- Test de normalité  : QQ-plot + Shapiro-Wilk ($H_0:$ Normalité).

```{r Tests statistique sur les residus}
# Test des hypotheses statistiques du modele :
# Autocorrelation, Durbin-Watson :
lmtest::dwtest(formula = modele_r, alternative = 'two.sided')
# Homoscedasticite, Breush-Pagan :
bptest(formula = modele_r)
# Normalite des residus, Shapiro : 
shapiro.test(modele_r$residuals)
```

### Résidus standardisés et stundentisés :

Soit $e_{i}$, le résidu pour la $i^{ème}$ observation : $e_{i}=y_{i}-\hat{y}_{i}$ 

Afin de rendre les résidus comparables, on divise par l'écart-type (*standardisation*). 

$$t_{i}=\frac{e_{i}}{\hat{\sigma}_{residus}^{2}}=\frac{e_{i}}{\hat{\sigma}\sqrt{(1-h_{i})}} \sim \mathcal{T}_{n-p-1}$$

Le résidu est généralement jugé suspect si : $t_{i}>2 \space ou \space 3$

Remarque :Ce type de résidus est aussi appelé **studentisé interne** 

ATTENTION : le résidu qui fait l'objet de l'analyse a servi à calculer les paramètres du modèle. La solution de ce problème réside dans le résidu studentisé (par opposition au résidu standardisé et aussi appelé **studentisé externe**) :

$$t_{i}^{\star}=\frac{y_{i}-\bar{y}_{(-i)}}{\hat{\sigma}_{(-i)}\sqrt{(1-h_{(-i)})}}$$

Un résidu studentisé est considéré trop grand typiquement s'il dépasse $2$ (fractile d'ordre $0.95$)


```{r Diagramme des resdisus standardises et studentises}
# Diagramme des residus studentises :
fit_res = data.frame(
  name = names(rstudent(modele_r)),
  fit = fitted(modele_r),
  rstud = rstudent(modele_r)
)
ggplot(fit_res, aes(x = fit, y = rstud)) +
  geom_point() + xlab("Fitted Values") +
  ylab("Studentized Residuals") +
  geom_hline(
    yintercept = 2,
    linetype = "dashed" ,
    color = "blue",
    linewidth = 0.8
  ) +
  geom_hline(
    yintercept = -2,
    linetype = "dashed" ,
    color = "blue",
    linewidth = 0.8
  ) +
  geom_text(
    data = subset(fit_res, abs(rstud) > 2),
    aes(x = fit, y = rstud, label = name),
    nudge_x = 3
  )
```

Retrouver les points déja observés à l'aide d'autres analyse oriente vers la conclusion que ces observations sont abérrantes. 

On doit passer par l'analyse des leviers mais **de toutes façons, observations ne pourront être retirées de l'échantillon qu'après une analyse métier**.

### Etude des leviers :

Le **levier** indique la distance au centre du nuage dans l'espace défini par les variables exogènes. C'est l'élément correspondant à l'observation sur la **diagonale principale de la Hat-matrix** :

Soit $h_{ii}=x_{i}(X^{T}X)^{-1}x_{i}^{T}$, le **levier** de l'observation $i$.

**Rappel :** La matrice $H$ (Hat-matrix) se présente comme suit :

$$
H =
\begin{pmatrix}
h_{11} & h_{12} & \cdots & h_{1n} \\
h_{21} & h_{22} & \cdots & h_{2n} \\
\vdots & \vdots & \ddots & \vdots \\
h_{n1} & h_{n2} & \cdots & h_{nn}
\end{pmatrix}_{n\times n}
$$

Dans la méthode des moindres carrés, les leviers mesurent l'influence d'une observation sur les estimations des coefficients de régression.

Il est généralement considéré haut pour l'observation $i$ si :

$$h_{ii} \geq 2\times \frac{P+1}{n}$$

```{r Etude des leviers}
# Detection de points aberrants / influents :
# Etude du levier :
ggplot(res, aes(seq_along(.hat), .hat)) + geom_col(width = 0.1, colour = "blue") +
  labs(x = "Observation", y = "Leverage") + geom_text(
    label = rownames(df),
    check_overlap = T,
    vjust = -0.8,
    size = 3
  ) +
  geom_hline(yintercept = 2 * length(coefficients(res)) / nrow(df),
             colour = "red")
# La liste des leviers est directement sur la diagonale de la matrice H :
# Il faut recalculer les matrices caracteristiques du model restreint :
X=as.matrix(cbind(1,modele_r$model))
H=X%*%solve(t(X)%*%X)%*%t(X)
# Puis extraire la diagonale de la matrice H :
leviers_modele_r = diag(H)
```

### Distance de Cook :

Mesure la distance entre le paramètre $\beta_{j}$ calculé avec l'observation $i$ et celle calculée sans cette observation : 

$$D_{i}=\frac{\sum_{i=1}^{n}(\hat{y}_{j}-\hat{y}_{j(-i)})}{\hat{\sigma}^{2}(P+1)}=\frac{h_{ii}\space e_{i}^{2}}{\hat{\sigma}^{2}(P+1)(1-h_{ii})^{2}}$$

Remarque : une distance de *Cook* suppérieur à $1$ est généralement considérée comme anormale, cependant d'autres valeurs de références sont possibles.


```{r mesures des leviers}
# Distance de Cook :temp_1=round(influence.measures(res)$infmat,2)
ols_plot_cooksd_chart(modele_r)
ols_plot_cooksd_bar(modele_r)
```

Pour l'interprêtation, on rapproche ces résultats de l'analyse des **résidus studentisés**. Les observations mises en évidence par cette méthode sont candidates à une analyse métier. 

## Autres mesures de l'influence : 

### Covratio :

Mesure la variabilité des estimations lorsqu'on retire une observation : 

$$Covratio = \frac{det\left(\hat{\sigma}^{2}\left(X_{(-i)}^{T}X_{(-i)}\right)^{-1}\right)}{det\left(\hat{\sigma}^{2}\left(X^{T}X\right)^{-1}\right)}$$

Propriétés : 

- $Covratio<1\space :$ l'observation dégrade l'estimation.

- $Covratio>1\space :$ l'observation améliore l'estimation.

- Limite généralement admise pour une observation influente : $\pm\frac{1+3(P+1)}{n}$

### DFBETAS :

Mesure standardisée de la différence entre l'estimation de $\beta_{j}$ avec et sans l'observation :

$$DFBETAS_{ij}=\frac{b_{j}-b_{j(-i)}}{\hat{\sigma}_{(-i)}^{2}\left(X^{T}X\right)_{(j+1),(j+1)}^{-1}}$$

Considérée importante pour : $|DFBETAS_{ij}|\geq \frac{2}{\sqrt{n}}$

### DFFITS :

Mesure standardisée de la distance entre les estimations $\hat{Y}$, avec et sans l'observation : 

$$DFFITS_{i}=\frac{\left(\hat{Y}_{i}-\hat{Y}_{i}(-i )\right)}{\hat{\sigma}_{(-i)}\sqrt{h_{ii}}}$$

Considérée importante pour : $|DFFITS_{i}|\geq 2\sqrt{\frac{(P+1)}{n}}$

```{r Autres mesures des influences}
# Autres mesures : 
Influences_modele_r=round(influence.measures(modele_r)$infmat,2)
# Labels des colonnes de temp_1 :
labels(Influences_modele_r)[[2]]
# Exemple de classement du resultat sur un colonne particuliere (ex. distance de cook) :
Influences_modele_r[order(as.data.frame(Influences_modele_r)$cook.d,decreasing = T),]
```

En rapprochant ces observations des autres faites précédement, on a une idée précise de l'influence de chaque observation sur les estimation des paramètres, des prédictions et de l'évolution de la variance.

## Références :

• Pierre-André Cornillon, Eric Matzner-Løber (2011). Régression avec R, Springer

• F. Housson (2018). R pour la statistique et la science des données (2ème ed.), PUR

• Ricco Rakotomalala. Pratique de la Régression Linéaire Multiple. Diagnostic et sélection de variables,
https://eric.univ-lyon2.fr/~ricco/cours/cours/La_regression_dans_la_pratique.pdf.
20
