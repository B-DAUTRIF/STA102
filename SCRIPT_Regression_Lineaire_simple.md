---
title: "RL_simple"
output: html_notebook
---

```{r Chargement des bibliotheques}
library(prettyR)
library(car)
library(ggplot2)
library(ggfortify)
library(lmtest)
library(olsrr)
```

## Données :

Les données sont chargées puis les variables pertinentes sont sélectionnées. L'affichage du nuage de points donne une première intuition de la relation entre $X$ et $Y$ :

```{r Chargement et mise en fore des donnees}
# Chargement du fichier de donnees :
df <- read.delim("D:/Etudes/2025_2026/STA102/Regression_simple/appartements.txt")

# Extraction des vecteurs de donnees :
y=df[2]
x=df[3]

# Concatenation des donnes pertinentes :
# La première colonne est un numero sans valeur statistique qui peut être # retiree : 
df=cbind(y,x)
names(df)=c("y","x")

# Visualisation du nuage de points :
plot(df)

# Visualisation du nuage de points avec GGPLOT :
# g <- ggplot(appartement, aes(x = Superficie , y = Prix)) + geom_point() + labs(x = "Superficie en m2", y = " Prix en milliers de Francs")
# 
# g

# Affichage d'informations statistiques generales sur le jeu de donnees :
describe(df)
```

### Corrélation : 

Coefficient de corrélation :

$$R=\frac{\sum_{i=1}^{n}(x_{i}-\bar{x})(y_{i}-\bar{y})}{\sqrt{\sum_{i=1}^{n}(x_{i}-\bar{x})\sum_{i=1}^{n}(y_{i}-\bar{y})}}$$

Matrice de corrélation : 

$$
R =
\begin{pmatrix}
1      & r_{12} & r_{13} & \cdots & r_{1n} \\
r_{21} & 1      & r_{23} & \cdots & r_{2n} \\
r_{31} & r_{32} & 1      & \cdots & r_{3n} \\
\vdots & \vdots & \vdots & \ddots & \vdots \\
r_{n1} & r_{n2} & r_{n3} & \cdots & 1
\end{pmatrix}
$$

Propriétés : 

- $r\in[-1,1]$

$$
r = \frac{s_{XY}}{s_X s_Y} \ \Longrightarrow \ 
s_{XY} = r \times s_X s_Y \ \Longrightarrow \
\frac{s_{XY}}{s_X} = r \times \frac{s_X s_Y}{s_X} \ \Longrightarrow \ 
\frac{s_{XY}}{s_X} = r \frac{s_Y}{s_X}
$$

**Attention :**

- $X \perp Y \Rightarrow r=0$

- $r=0 \ \nRightarrow X \perp Y$

```{r Correlations}
# Formule du calcul du coefficient de correlation r :
# sqrt(cov(y,x)/(var(x)*var(x)))

# Representation graphique des correlations entre variables :
pairs(df)
# Matrice de correlations :
cor(df)
# Matrice graphique de correlations :
corrplot::corrplot(cor(df))
```

## Modèle : 

Compte tenu de la forme du nuage et données dont on dispose, on décide si on peut appliquer un modèle du type :

$$Y_{i}=\alpha+\beta x_{i}+\epsilon_{i}$$

On créé donc le modèle et on affiche la droite de régression :

#### Equation de la droite des moindres carrés :

$$\hat{y}=a+bx=\bar{y}-b\bar{x}+r\frac{S_{Y}}{S_{X}}x=\bar{y}-r\frac{S_{Y}}{S_{X}}\bar{x}+r\frac{S_{Y}}{S_{X}}x=\bar{y}+r\frac{S_{Y}}{S_{X}}(x-\bar{x})$$

Propriétés : 

- La droite des moindres carrés passe par $(\bar{x} \ ; \ \bar{y})$, qui est le centre de gravité du nuage.       
- Le signe de sa pente est le même que celui de $S_{XY}$.

**Rappel :** $\forall \space i \space : \space \epsilon_{i}\sim \mathcal{N}{(0 \space ;\space \sigma^{2})}$

#### Propriétés de $A$ et $B$ :

- **Linéarité** car ils sont combinaisons linéaires des $Y_{i}$ :

$$\delta_{i}=\sum^{n}_{i=1}\frac{(x_{i}-\bar{x})}{\sum^{n}_{i=1}(x_{i}-\bar{x})^{2}}Y_{i.} \ \qquad \Longrightarrow \qquad B=\Sigma^{n}_{i=1}\delta_{i}Y_{i} \qquad et \qquad A=\Sigma^{n}_{i=1}\left(\frac{1}{n}\bar{x}\delta_{i}\right)Y_{i}$$

- Absence de **biais** : 

$$E[\ B \]=E\left[\frac{\sum_{i=1}^{n} (Y_i-\bar Y)(x_i-\bar x)}{\sum_{i=1}^{n}(x_i-\bar x)^{2}}\right]=  \frac{\sum^{n}_{i=1} E[Y_{i}-\bar Y](x_{i}-\bar x )}{\sum_{i=1}^{n}(x_{i}-\bar{x})^{2}}$$

On sait que :

$$E[\ Y_{i} \]=E[\ \alpha+\beta x_{i}+\epsilon_{i} \]=\alpha+\beta x_{i}$$

$$E[\ \bar{Y} \ ]=\frac{1}{n}\Sigma_{i}E[\ Y_{i} \ ]=\frac{1}{n}\Sigma_{i}(\alpha+\beta x_{i})=\alpha+\beta \bar{x}$$

$$\Rightarrow E[\ B \ ]=  E\left[\ \frac{\sum^{n}_{i=1} \beta (x_{i}-\bar{x})(x_{i}-\bar{x})}{\sum^{n}_{i=1}(x_{i}-\bar{x})^{2}} \ \right]=\beta$$

$$
E[\ A\ ]=E[\ \bar{Y}-\beta\bar{x}\ ]=E[\ Y\ ]-\beta x=(\alpha+\beta x+E[\ \epsilon\ ])-\beta x=\alpha
$$

On démontre que : 

$$
V[\ B\ ]=\frac{\sigma^{2}}{nS^{2}_{X}}
$$

$$
V[\ A\ ]=\frac{\sigma^{2}}{n}\left(1+\frac{\bar{x}}{S^{2}_{X}}\right)
$$

Où $\sigma^{2}$ est la variance des résidus et $S_{X}^{2}$ celle de $X$ : $S_{X}=\frac{1}{n}\Sigma_{i}(x_{i}-\bar{x})^{2}$

- **non indépendance** : $A \not\perp B$

On démontre que : 

$$Cov(A \ ; \ B)=\frac{-\bar{x}\sigma^{2}}{\Sigma^{n}_{i=1}(x_{i}-\bar{x})^{2}}$$

$$\rho(A \ ; \ B)=\frac{-\bar{x}}{\sqrt{\Sigma^{n}_{i=1}(x_{i}^{2})/n}}$$

- $A$ et $B$ sont des estimateurs de **variance minimale** (théorème de *Gauss-Markov*).

```{r Construction du modele}
# On construit le modele lineaire de regression simple :
modele=lm(y~x,data=df)

# Visualisation de la droite de regression :
plot(y~x,data=df)
abline(modele, col='red')
```

### Paramètres du modèle : 

$$a=\bar{y}-b\bar{x}$$

$$
b= \frac{\sum_{i=1}^n (x_i-\bar x)(y_i-\bar y)}{\sum_{i=1}^n (x_i-\bar x)^2}
$$

```{r Rapports sur le modele}
# Synthese :
summary(modele)
# Coefficients avec IC :
Confint(modele)
```

## Equation de la droite de régression : 

$$\hat{Y}_{i}=-147.33+26.77x_{i}$$

```{r Affichage avance}
# Fonction d'affichage de la droite de régression utilisant GGPLOT :
ggplot(df, aes(x = x, y = y)) + geom_point() + 
  geom_smooth(method = lm, se = F) +
  geom_text(
    label = row.names(df),
    vjust = -1,
    check_overlap = TRUE,
    size = 3
  )
```

Remarques sur les résultats de l'exemple : 

- L'intercept n'a pas ici de valeur en soi puisqu'il n'est pas envisageable de disposer d'un appartement de $0m^{2}$.          
- L'écart aux estimations va croissant sur $X$, ce qui laisse supposer un éventuel problème d'**hétéroscédasticité**.          
- Certaines observations ($18,25,6,4$) montrent un résidu important. On doit supposer (et vérifier) des observations abérrantes.          
- Une observation ($16$) montre une valeur de $X$ particulièrement grande. La question se posera de sa pertinence ou d'une observation abérrante.          

## Table d'analyse de la variance : 

```{r Table analyse de variance}
# On affiche les differentes sommes de carres du modele :
anova(modele)
```

## Intervalles de confiance (IC) :

On vérifie les doutes sur les observations en observant les intervalles de confiance :

```{r Intervalles de Confiance}
# Creation d'un jeu de donnees incluant les predictions et les
# intervalles de confiance :
IC = as.data.frame(cbind(df, predict(
  modele, interval = "confidence", level = 0.95
)))
```

Puis on affiche les IC sur la droite de régression :

```{r Droite de regression avec IC}
# Affichage de la droite de regression avec IC :
ggplot(appartement, aes(x=Superficie, y=Prix))+ geom_point()+geom_smooth(method=lm, se=T)
```

Remarque : On remarque un élargissement des IC a mesure que les valeurs s'écartent de $\bar{X}$.

## Prédictions :

Etude du cas où l'on ne connait que la variable $X$ à partir de laquelle on souhaite **inférer** $Y$ :

Ci-dessous un exemple pour $X=150$

```{r Exemple Intervalle de Prediction}
predict(modele,list(x=150),interval = "prediction")
```

On peut étendre cette démarche au jeu de données entier (cf. le logiciel nous prévient que nous allons utiliser les mêmes données que celles qui ont servi à estimer le modèle) :

```{r Intervalle de Prediction sur le jeu de donnes}
# Estimation des predictions sur l'ensemble du jeu de donnes :
IP = as.data.frame(cbind(
  Prix = df$y,
  Superficie = df$x,
  predict(modele, interval = "prediction")
))
```

On affiche l'intervalle de prédiction sur la droite de régression : 

```{r Itervalles de Prediction sur la droite de regression}
# Affichage de l'intervalle de prediction sur la droite de regression :
ggplot(IP, aes(x = Superficie, y = Prix)) +
  geom_line(aes(y = lwr), color = "red", linetype = "dashed") +
  geom_line(aes(y = upr), color = "red", linetype = "dashed") +
  geom_point() + 
  geom_smooth(method = lm, se = T)+ 
  geom_text(label = row.names(df),
vjust = - 1, check_overlap = TRUE, size = 3)
```

## Etude des résidus :

On affiche d'abord la synthèse détaillée par observations du modèle  (variables, estimations, résidus, leviers, variances, distance de Cook) :

```{r synthese detaillee du modele}
# Synthèse detaillee (levier, fitted, hat, cook) :
broom::augment(modele)
```

On visualise le tracé des résidus contre les prédictions :

```{r Graphique des residus et valeurs estimees}
# Partitionnement de l'ecran pour l'affichage :
split.screen(1:2)

# On peut visualiser les valeurs evaluees du modele :
screen(2) ; plot(modele$fitted.values)
# On peut visualiser les residus du modele :
screen(1) ; plot(modele$residuals)

close.screen(all = TRUE)
```

```{r Graphique des residus VS valeurs estimees}
# Utilisation de la bibliothéque GGFORTIFY :
autoplot(
  modele,
  which = 1,
  ncol = 1,
  label.size = 3,
  label.hjust = -0.8,
  label.n = 6
)
```

### Evaluation des hypotheses sur les residus :

**Rappel :** 

- Indépendance : $\forall \space i\neq j \space : \space \epsilon_{i}\perp \epsilon_{j}$

- Homoscédasticité et identité de loi : $\forall \space i\space : \space \epsilon_{i}\sim \mathcal{N}(0 \space ; \space \sigma^{2})$

### Recherche d'autocorrélation : 

#### Test de Durbin-Watson :

$H_{o} : coefficient =0$

Fonction implémentée dans la librairie `DescTools`

```{r Test autocorrelation}
# Test de Durbin-Watson : 
durbinWatsonTest(modele,alternative = "two.sided")

```

### Graphique standards des résidus du modèle :

4 graphiques de base sont disponibles concernant les résidus du modèle : 

- Résidus ~ esimations

- Quantile-Quantile des résidus : en abscisse les résidus standardisés, en ordonnée les quantiles théoriques issus du modèle.

- Estimations ~ résidus standardisés

- Résidus ~ levier

```{r Graphiques standards des residus du modele}
# Lot de graphiques standards des residus du modele :
plot(modele)
```

Voire la librairie **OLSRR** pour une liste plus complète d'indicateurs et les graphiques correspondants.

### Graphique des valeurs absolues des résidus studentisés internes :

Recherche l'augmentation de valeurs des résidus.

```{r Graphique residus studentises internes VA}
# Fonction graphique utilisant GGPLOT pour les valeurs absolues 
# des résidus studentisés internes :
autoplot(modele, which = 3, ncol = 1, label.size = 3, label.hjust = -0.8, label.n = 6)
```

### Graphique des résidus studentisés internes :

Elimination de l'échelle pour disposer de valeurs de référence (typiquement $2, \space 2.5$ ou $3$)

```{r Graphique residus studentises internes}
# Fonction graphique utilisant GGPLOT pour les residus studentises
# internes :
autoplot(
  modele,
  which = 3,
  ncol = 1,
  label.size = 3,
  label.hjust = -0.8,
  label.n = 6
) +
  aes(.fitted, .stdresid) + labs(x = "fitted values", y = "standardized residuals") +
  ggtitle("")
```

Test d'héthérogénéïté de Breusch-Pagan :
(voir wiki : Régression linéaire simple)

```{r Test heterogeneite de Breusch Pagan}
# Test d'hetherogeneite de Breusch-Pagan :
bptest(modele)
```

### Evaluation de la normalité des résidus :

Affichage du graphique des rapports de quantiles (QQplot) :

```{r Graphique Quatile VS Quantile}
# Quantile-Quantile-plot (QQ-plot) :
autoplot(
  modele,
  which = 2,
  ncol = 1,
  label.size = 3,
  label.hjust = -0.8,
  label.n = 6
)
```

### Test de normalité de Shapiro-Wilk :
```{r Test Shapiro Wilk}
shapiro.test(residuals(modele))
```










