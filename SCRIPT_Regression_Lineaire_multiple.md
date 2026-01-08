---
title: 'Régression linéaire multiple'
output:
  html_document:
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

res = modele_r

# Diagramme des residus studentises :
fit_res = data.frame(
  name = names(rstudent(res)),
  fit = fitted(res),
  rstud = rstudent(res)
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

# Diagramme des leviers :
ggplot(res, aes(seq_along(.hat), .hat)) + geom_col(width = 0.1, colour = "blue") +
  labs(x = "Observation", y = "Leverage") + geom_text(
    label = rownames(df),
    check_overlap = T,
    vjust = -0.8,
    size = 3
  ) +
  geom_hline(yintercept = 2 * length(coefficients(res)) / nrow(df),
             colour = "red")
```

Retrouver les points déja observés à l'aide d'autres analyse oriente vers la conclusion que ces observations sont abérrantes. 

On doit passer par l'analyse des leviers mais **de toutes façons, observations ne pourront être retirées de l'échantillon qu'après une analyse métier**.

```{r }
# Detection de points aberrants / influents :
# Etude du levier :
ggplot(res, aes(seq_along(.hat), .hat)) + geom_col(width = 0.1, colour = "blue") +
  labs(x = "Observation", y = "Leverage") + geom_text(
    label = rownames(ozone),
    check_overlap = T,
    vjust = -0.8,
    size = 3
  ) +
  geom_hline(yintercept = 2 * length(coefficients(res)) / nrow(ozone),
             colour = "red")

# Distance de Cook :temp_1=round(influence.measures(res)$infmat,2)
ols_plot_cooksd_chart(res)
ols_plot_cooksd_bar(res)

# Autres mesures : 
temp_1=round(influence.measures(res)$infmat,2)
# Labels des colonnes de temp_1 :
labels(temp_1)[[2]]
# Exemple de classement du resultat sur un colonne particuliere (ex. distance de cook) :
temp_1[order(as.data.frame(temp_1)$cook.d,decreasing = T),]
```

$$
X_{n \times p},
\quad
(X^\top X)^{-1}_{p \times p},
\quad
H_{n \times n} = X_{n \times p} (X^\top X)^{-1}_{p \times p} X^\top_{p \times n}
$$

$$
H = (h_{ij})_{1 \le i,j \le n},
\qquad
h_{ij} = x_i^\top (X^\top X)^{-1} x_j
$$
$$
H = X (X^\top X)^{-1} X^\top
$$

$$
H \in \mathbb{R}^{n \times n}, \quad
X \in \mathbb{R}^{n \times p}, \quad
X^\top X \in \mathbb{R}^{p \times p}
$$
