---
title: "Régression linéaire multiple"
output: html_notebook
---

Ce bloc-note reprend le contenu du cours STA102 du Conservatoire National des Arts et Métier (CNAM) sur les **modéles linéaires**. Son but est de décrire les étapes d'une **régression linéaire multiple** à *P* variables explicatives, en rassemblant à la fois le matériel théorique et pratique pour travailler efficacement sur des données. 

Ce script est destiné à un usage personnel. J'essaie au maximum qu'il puisse être utilisé tel quel pour effectuer les calculs et la production des différents résultats (matrices, vecteurs paramètres et de variance, testes statistiques) : a condition pour cela d'effectuer les changement nécessaires pour l'import de nouveaux jeux de données. 

Un éventuel lecteur intéressé par ce document est invité à laisser suggestions et remarques (par le biais de GitHub).

### Préparation des données et de l'environnement : 

Chargement des libraries :
```{r Chargement des bibliotheques}
library('corrplot')
library('prettyR')
```

Chargement du jeu de donnée :
```{r Chargement du jeu de donnees}
df <- read.csv2(
  "D:/Etudes/2025_2026/STA102/3.Regression multiple/TP_regression_multiple/ozone.csv"
)
```

description du jeu de données :
```{r Description du jeu de donnees}
# Description des variables : 
# summary(df)
# sum(is.na(df))
# dim(df)
# sum(is.na(df))
describe(df)
```

### Création des matrices du modèle :

Le **modèle général** de la régression linéaire peut s'écrire de la façon suivante : 

$$y_{i}=\beta_{0}+\Sigma_{j=1}^{p} \ \beta_{j} \ x_{ij} + \epsilon_{i}$$

```{r Preparation des matrices}
# Vecteur Y (variable a expliquer) :
Y = as.matrix(df[, 1])
# Matrice X (intercept + observations) :
X = cbind(matrix(
  nrow = dim(df)[1],
  ncol = 1,
  data = 1
), as.matrix(df[, -1]))
# Vecteur vide B (paramètres du modele) : 
B = matrix(nrow = dim(df)[2], ncol = 1)
# Vecteur vide Epsilon (erreurs aleatoires) :
Epsilon = matrix(nrow = dim(df)[1], ncol = 1)
```

Vérification de la taille des matrices :
```{r Banc de verification des tailles de matrices}
# Verification des dimensions des matrices : 
# print("Dimensions de Y :")
# dim(Y)
# print("Dimensions de X :")
# dim(X)
# print("Dimensions de B :")
# dim(B)
# print("Dimensions de Epsilon :")
# dim(Epsilon)
```

Affichage de la matrice de corrélation :
```{r matrice de correlation des observations}
# Affichage de la matrice de correlation de X :
corrplot(cor(X[,-1]))
```

### Ecriture matricielle :

Le modèle peut maintenant s'écrire de la façon suivante :

$$\mathcal{y}=X\beta +\epsilon$$

Propriétés :

$$\epsilon \sim \mathcal{N}(0,\sigma^{2}I)$$      

$\Rightarrow \ y \sim \mathcal{N}(X\beta,\sigma^{2}I)$      
$\Rightarrow \ (y-X\beta) \sim \mathcal{N}(0,\sigma^{2}I)$      
$\Rightarrow \ b \sim \mathcal{N}(\beta,\sigma^{2}(X^{T}X)^{-1})$      

### Estimation des paramètres du modèle :

Méthode des moindres carrés :

$$Min((y-X\beta)^{T}(y-X\beta))$$

Alors :

$$(y-X\beta)^{T}(y-X\beta)=(y^{T}-\beta^{T}X^{T})(y-X\beta)=y^{T}y-\beta^{T}X^{T}y-y^{T}X\beta+\beta^{T}X^{T}X\beta$$      

$(\beta^{T}X^{T}y=y^{T}X\beta)$ est un **scalaire**, on peut donc écrire : 

$$y^{T}y-\beta^{T}X^{T}y-y^{T}X\beta+\beta^{T}X^{T}X\beta=y^{T}y-2\beta^{T}X^{T}y+\beta^{T}X^{T}X\beta$$

Dérivation de $(y-X\beta)^{T}(y-X\beta)$ par $\beta$:      

$$\frac{\partial(y^{T}y-2\beta^{T}X^{T}y+\beta^{T}X^{T}X\beta)}{\partial\beta}=2X^{T}y-2X^{T}X\beta=2X^{T}(y-X\beta)$$

Résolution :      

$$\Rightarrow 2X^{T}y-2X^{T}X\beta=0$$ 
$$X^{T}y-X^{T}X\beta=0$$
$$X^{T}X\beta=X^{T}y$$
$$\beta=(X^{T}X)^{-1}X^{T}y$$

### Estimateur des moindres carrés (EMC) :

$$EMC=A=(X^{T}X)^{-1} \ X^{T}$$ 
$$b=(X^{T}X)^{-1} \ X^{T} \ Y=AY$$

**Attention :** cette solution exige que la matrice $X$ soit inversible, il est donc **obligatoire que les variables explicatives soient indépendantes**. Soit :   

- aucune variable explicative ne doit être combinaison linéaire des autres      
- X doit être de rang plein

```{r Estimateur des moindres carres}
# Estimateur des moindres carres (EMC) :
EMC <- function (Y, X, A_ret = F) {
  # Entree : vecteur Y, matrice X, option de retour de l'EMC par la fonction
  # - calcul l'EMC (A) a partir de la matrice X
  # - applique l'estimateur a la matrice d'observations X
  # - selon l'option choisie, ajoute au l'EMC au retour de la fonction
  # Sortie : liste comprenant le vecteur des estimations (Y_hat) +/- l'EMC
  A = solve(t(X) %*% X) %*% t(X)
  B = A %*% Y
  if (A_ret == T) {
    res = list(B, A)
  }else{
    res = B
  }
  return(res)
}

# Realisation de l'EMC :
B = EMC(Y = Y, X = X, A_ret = F)
# Estimation de Epsilon :
Epsilon = Y - X %*% B
```

### Propriétés :      
- **linéarité** : $A$ étant combinaison linéaire des données observées
- **Sans bias** :

$$b=(X^{T}X)^{-1} \ X^{T} \ Y=(X^{T}X)^{-1} \ X^{T}(X\beta+\epsilon)=(X^{T}X)^{-1} \ X^{T}X\beta+(X^{T}X)^{-1}X^{T}\epsilon=\beta+A\epsilon$$

Avec :

$$E[ \ b \ ]=E[ \ \beta+A\epsilon \ ]=\beta+AE[ \ \epsilon \ ]=\beta$$

Et $V[ \ b \ ]$, la matrice de covariance de $b$ :

$$V[ \ b \ ]=V[ \ Ay \ ]=A^{T} \ V[ \ y \ ] \ A=\sigma^{2}A \ A^{T}$$
$$=\sigma^{2}(X^{T}X)^{-1}X^{T}X(X^{T}X)^{-1}=\sigma^{2}(X^{T}X)^{-1}$$

Parmi tous les estimateurs lin´eaires non biais´es de β, b est l’estimateur de variance minimale (*théorème de Gauss-Markov*).

## Matrices caractéristiques :

- **Hat-matrice (H) :**

$$\hat{y}=X\beta=X(X^{T}X)^{-1}X^{T}y=Hy$$
$$H=X(X^{T}X)^{-1}X^{T}$$

Propriétés :       
- matrice de projection dans $\mathbb{R^n}$ sur le sous-espace engendré par les vecteurs représentant les variables explicatives.      
- $H_{n\times n}$      
- $H^{T}=H$      
- $H^{2}=H$      
- $Rang(H)=Rang(X)=P+1$      
- $h_{ii}=x_{i}(X^{T}X)^{-1}x_{i}^{T}$, est le **levier** de l'observation $i$.

```{r Hat matrice}
matrice_Hat <- function(X) {
  H = X %*% solve(t(X) %*% X) %*% t(X)
  return(H)
}
H = matrice_Hat(X)
leviers = diag(H)
```

- **matrice M :**

$$e=y-\hat{y}=y-Hy=(I-H)y$$
$$M=(I-H)=I-(X^{T}X)^{-1}X^{T}$$      

On note : $e=My$

Propriétés :      
- $M^{T}=M$      
- $M^{2}=M$      
- $M$ est une matrice de projection orthogonale.      
- $M$ est orthogonale à $H$.      
- $Rang(M)=Trace(M)=n-p-1$      

```{r M matrice}
matrice_M <- function(X) {
  H = matrice_Hat(X)
  M = diag(dim(H)[1]) - H
  return(M)
}

M = matrice_M(X)
# Verification de M par comparaison des residus observes (a 6 decimales) :
residus_obs = M %*% Y
if(all(round(residus_obs,6)==round(Epsilon,6))){
  print("Matrice M correcte")
}

V_Y = as.numeric(var(Y))
V_e = V_Y * H
```

### Espérance et Variance de $(y-X\beta)$

$$e=My=M(X\beta+\epsilon)=(I-X(X^{T}X)^{-1}X^{T})X\beta+M\epsilon=M\epsilon$$
$$\Rightarrow E[\ y-X\beta\ ]=E[\ M\epsilon\ ]=ME[\ \epsilon\ ]=0$$

et : 

$$V[\ X\beta+\epsilon \ ]=V[\ M\epsilon \ ]=M^{2}V[\ \epsilon \ ]=\sigma^{2}M$$

```{r Loi des residus observes}
# Verification de la loi des residus observes :
# Moyenne (doit etre nulle)
mean(Epsilon)
# Estimation de la variable des residux observes : 
variance_residus <- function(X, Y, Epsilon_ret = F) {
  M = matrice_M(X)
  residus_obs = M %*% Y
  # On choisit ici l'estimateur non biaise de sigma_2 :
  sigma_2 = sum(residus_obs ** 2) / (dim(X)[1] - dim(X)[2] - 1)
  # option de renvoi du vecteur des residus :
  if (Epsilon_ret == T) {
    res = list(sigma_2, residus_obs)
  } else{
    res = sigma_2
  }
  return(res)
}
sigma_2 = variance_residus(X, Y)

# Estimation de la variance des estimateurs de parametres (B) :
variance_betas <- function(X, Y, Cov_ret=F) {
  # Preparation d'une matrice de resultats :
  res = matrix(nrow = dim(X)[2])
  # Extraction du vecteur de residus :
  residus_obs = M %*% Y
  # Calcul de la variance des residus (sigma_2) :
  variance_residus = sum(residus_obs ** 2) / (dim(X)[1] - dim(X)[2] - 1)
  # Calcul de la matrice de variance-covariance des estimateurs :
  V_B = variance_residus * solve(t(X) %*% X)
  # Selection de la variance des estimateurs Betas :
  for (i in 1:dim(X)[2]) {
    res[i] = V_B[i, i]
  }
  if (Cov_ret==T){
    res=list(res,V_B)    
  }
  return(res)
}
V_Betas=variance_betas(X,Y)
```

#### Estimation par Maximum de de Vraissemblance (EMV) :

**Rappels :** 

- $\epsilon \sim \mathcal{N}(0,\sigma^{2}I)$

- Pour l'**espérance** l'$EMV$ coïncie à l'$EMC$.

- Pour la **variance** $\sigma^{2}$ :

$$S^{2}_{n}=\frac{1}{n}\Sigma^{n}_{i=1}(y_{i}-\hat{y}_{i})^{2}=\frac{SCR}{n}$$

Avec :

$$E[\ SCR\ ]=E[\  (y-X\beta)^{T}(y-X\beta)\ ]=E[\  (Me)^{T}(Me)\ ]$$
$$=E[\  e^{T}M^{T}Me\ ]=E[\  e^{T}Me\ ]=tr[\  e^{T}Me\ ]=tr[\  M\ ] \ E[\  e^{T}e\ ]$$
$$=\sigma^{2} \ tr[\ M\ ]=\sigma^{2}(n-p-1)$$

L'estimateur **non biaisé** de $\sigma^{2}$ est donc : 

$$S^{2}_{n-p-1}=\frac{(y-X\beta)^{T}(y-X\beta)}{n-p-1}=\frac{1}{n-p-1}\Sigma_{i=1}^{n}(y_{i}-B_{0}-B_{1}x_{i1}...-B_{p}x_{ip})^{2}$$
$$\Rightarrow E[\ S^{2}_{n-p-1} \ ]=E \left[\frac{(y-X\beta)^{T}(y-X\beta)}{n-p-1} \ \right]=\frac{\sigma^{2}(n-p-1)}{n-p-1}=\sigma^{2}$$

**Rappel :** $V[\ b \ ]=\sigma^{2}(X^{T}X)^{-1}$ 

Donc : 

$$\hat{\sigma}^{2}_{b}=\hat{\sigma}^{2}(X^{T}X)^{-1}$$
$$\Rightarrow \ \hat{\sigma}^{2}_{b_{j}}=\hat{\sigma}^{2}[\ (X^{T}X)^{-1}\ ]_{(j+1,j+1)}$$

**Rappel :** $b=Ay$ 

$B$ est donc lui aussi un vecteur gaussien, comme transformée linéaire d'un autre vecteur gaussien. On peut lui appliquer des tests statistiques :

$$\frac{B_{j}-\beta_{j}}{\hat{\sigma}[\ (X^{T}X)^{-1/2}\ ]_{(j+1,j+1)}}\sim \mathcal{N}(0,1) \ ; \ j=0,...,p$$

### La régression comme projection orthogonale :

*D'un point de vue géométrique, la régression est une projection orthogonale de Y dans l'espace engendré par les variables explicatives.*

Soit le triangle $\bar{y} y \hat{y}$, rectangle en $\hat{y}$ :

$$SCT=\Sigma_{i=1}^{n}(y_{i}-\bar{y})^{2}=||y-\bar{y}||^{2}$$
$$SCR=\Sigma_{i=1}^{n}(y_{i}-\hat{y}_{i})^{2}=||y-X\beta||^{2}=||y-\hat{y}||^{2}$$
$$SCE=\Sigma_{i=1}^{n}(\hat{y}_{i}-\bar{y})^{2}=||\hat{y}-\bar{y}||^{2}$$

On retrouve donc, grâce au **théorème de Pythagore**, la décomposition de la somme des carrés totale, entre somme des carrés expliqués et résiduels :

$$||y-\bar{y}||^{2}=||y-\hat{y}||^{2}+||\hat{y}-\bar{y}||^{2} \ \Longrightarrow \ SCT=SCR+SCE$$

Lien avec le coefficient de **corrélation linéaire**. 

Soit $\eta$ l'angle formé par $(y- \bar{y})$ et $(\hat{y}-\bar{y})$ :

$$R^{2}=\frac{SCE}{SCT}=\frac{||\hat{y}-\bar{y}||^{2}}{||y-\bar{y}||^{2}}=Cos^{2}(\eta)$$
$$R=Cos(\eta)=\frac{||\hat{y}-\bar{y}||}{||y-\bar{y}||}$$

### Lois des sommes de carrés :

Soit le triangle ($X\beta,\hat{y},y$) :

$$y=X\beta + \epsilon \ \ \Rightarrow \ \ ||y-X\beta||^{2}=||\epsilon||^{2}=\epsilon^{T}\epsilon$$

$$Xb + \epsilon=M\epsilon \ \ \Rightarrow \ \ (y-Xb)^{T}(y-Xb)=(M\epsilon)^{T}(M\epsilon)=\epsilon^{T}M^{T}M\epsilon=\epsilon^{T}M\epsilon$$
$$||Xb + X\beta||^{2} = ||y-M\epsilon -(y-\epsilon)||^{2}=||(I-M)\epsilon||^{2}=||H\epsilon||^{2}=\epsilon^{T}H\epsilon$$

D'après le théorème de **Pythagore** :

$$||y-X\beta||^{2}=||y-Xb||^{2}+||Xb-X\beta||^{2}$$
$$\Rightarrow \ \epsilon^{T}\epsilon=\epsilon^{T}M\epsilon+\epsilon^{T}H\epsilon$$

Et comme :       
- $\epsilon \sim \mathcal{N}(0\ ; \ \sigma^{2}I) \ \Rightarrow \ \epsilon^{T}\epsilon\sim \sigma^{2}\chi^{2}_{n}$      
- $H$ et $M$ sont des matrices symétriques      
- $Rang(H)+Rang(M)=p+1+n-p-1=n$      

Alors : 

$$SCR\sim\chi^{2}_{(n-p-1)}$$

### Lois des côtés du triangle $(y \hat{y} \bar{y})$

On s'intéresse d'abord au triangle $(\bar{y},0,y)$, rectangle en $\bar{y}$ :

$$||y||^{2}=||\bar{y}||^{2}+||\bar{y}-y||^{2}$$
$$\Rightarrow \ ||y||^{2}=||\bar{y}||^{2}+SCT$$

Soit l'**opérateur de projection orthogonale** sur le vecteur $1$, tel que : 

$$\bar{y}=\frac{11^{T}}{n}y$$

Propriétés : 

- $Rang\left(\ \frac{11^{T}}{n} \right)=trace\left(\ \frac{11^{T}}{n} \right)=\Sigma^{n}_{i=1}\frac{1}{n}=1$

Donc : 

$$||\bar{y}||^{2}=y^{T}\left(\frac{11^{T}}{n}\right)\left(\frac{11^{T}}{n}\right)y=y^{T}\left(\frac{11^{T}}{n}\right)y$$

Et :

$$y-\bar{y}=y^{T}-\left(\frac{11^{T}}{n}\right)^{T}\left(\frac{11^{T}}{n}\right)y=\left(I-\frac{11^{T}}{n}\right)y$$
$$\Rightarrow \ ||y-\bar{y}||^{2}=y^{T}\left(I-\frac{11^{T}}{n}\right)^{T}\left(I-\frac{11^{T}}{n}\right)y=y^{T}\left(I-\frac{11^{T}}{n}\right)y$$

**Rappel :** $\hat{y}=Hy$

$$SCE=||\hat{y}-\bar{y}||^{2}=||Hy-\frac{11^{T}}{n}y||^{2}=||\left(H-\frac{11^{T}}{n}\right)y||^{2}=y^{T}\left(H-\frac{11^{T}}{n}\right)y$$

Donc :

$$||y||^{2}=||\bar{y}||^{2}+||y-\bar{y}||^{2}=y^{T}\left(\frac{11^{T}}{n}\right)Y+y^{T}\left(I-\frac{11^{T}}{n}\right)Y$$

Enfin :

$$SCT=SCR+SCE \ \Rightarrow \ y^{T}\left(\frac{11^{T}}{n}\right)y=y^{T}\left(I-H\right)y+y^{T}\left(H-\frac{11^{T}}{n}\right)y$$

Propriétés :

- $Rang(I)=n$      
- $Rang\left(\frac{11^{T}}{n}\right)=1$      
- $Rang\left(I-\frac{11^{T}}{n}\right)=n-1$      
- $Rang\left(I-H\right)=n-p$      
- $Rang\left(H-\frac{11^{T}}{n}\right)=p$      

**Rappels :**

- $\{y_{1},...,y_{n}\}$ est une suite de variables normales $iid$      
- $SCR \sim \chi^{2}_{(n-p-1)}$      

En appliquant le théorème de Cochrane au triangle $(y\hat{y}\bar{y})$, on obtient alors que : 

$SSi \ \ \{\forall i \in p : \beta_{i}=0\}$

$$SCT=\Sigma^{n}_{i=1} (y_{i}-\bar{y})^{2} \sim \sigma^{2} \ \chi^{2}_{n-1}$$

$$SCE=y^{T}\left(H-\frac{11^{T}}{n}\right)y \sim \sigma^{2}\chi^{2}_{p}$$

Et :

$$SCE \perp SCR$$

Cela justifie le teste de significativité globale du modèle.

## Tests d'hypothèses :

### Significativité globale :

$H_{0} : y_{i}=\beta_{0}+\epsilon_{i}$   

$H_{(1)} : y_{i=}\beta_{0}+\beta_{1}x_{i1}+...+\beta_{j}x_{ij}+...+\beta_{p}x_{ip}+\epsilon_{i}$

**Statistique :** 

$$F=\frac{SCE/p}{SCR/(n-p-1)} \sim \mathcal{F_{(p \ ; \ n-p-1)}}$$

*Rappel : sous* $H_{0}$ : 

- $SCE \sim \chi^{2}_{p}$
- $SCR \sim \chi^{2}_{n-p-1}$

Test associé à la **table d'analyse de la variance** :

| Source   | DDL         | Somme des carrés                                                                 | Carrés moyens                     | Statistique du test                                         | P-value                     |
|----------|-------------|----------------------------------------------------------------------------------|-----------------------------------|--------------------------------------------------------------|-----------------------------|
| Modèle   | $p$         | $\mathrm{SCE} = \sum_{i=1}^n (\hat y_i - \bar y)^2$                               | $\dfrac{\mathrm{SCE}}{p}$         | $F = \dfrac{\mathrm{SCE}/p}{\mathrm{SCR}/(n-p-1)}$           | $\mathbb{P}(F > F_{\text{obs}})$ |
| Résidu   | $n-p-1$     | $\mathrm{SCR} = \sum_{i=1}^n (y_i - \hat y_i)^2$                                 | $\dfrac{\mathrm{SCR}}{n-p-1}$     |                                                              |                             |
| Total    | $n-1$       | $\mathrm{SCT} = \sum_{i=1}^n (y_i - \bar y)^2$                                   | $\dfrac{\mathrm{SCT}}{n-1}$       |                                                              |                             |

### Modèles imbriqués :

Soit le modèle complet à $p$ paramètres et un modèle imbirqué à $r$ paramètres :

$H_{0}:y_{i}=\beta_{0}+\beta_{1}x_{i1}+...+\beta_{j}x_{ij}+...+\beta_{r}x_{ir}+\epsilon_{i}$      
$H_{1}:y_{i}=\beta_{0}+\beta_{1}x_{i1}+...+\beta_{j}x_{ij}+...+\beta_{p}x_{ip}+\epsilon_{i}$

*Remarque : sous* $h_{0} \ \rightarrow \ \{\forall j\in[r;p] \ : \ \beta_{j}=0\}$

Statistique : 

$$F=\frac{([SCE \ | \ H_{1}]-[SCE \ | \ H_{0}])/(p-r)}{[SCR \ | \ H_{1}]/(n-p-1)} \sim \mathcal{F}_{(p-r \ ;\ n-p-1)}$$

### Test de significativité des paramètres : 

$H_{0}: \ \beta_{j}=0$      
$H_{0}: \ \beta_{j} \neq 0$

Statistique : 

$$T=\frac{B_{j}}{S_{B_{j}}} \sim \mathcal{T}_{n-p-1}$$

## Critères de sélection des modèles : 

- Test de Fisher : donne la significativité de la différence entre modèles       
- $R^{2}_{a}$ : à la différence du $R^{2}$, permet de comparer des modèles dont le nombre de paramètres sont différents :

$$R^{2}_{a}=\frac{R^{2}(n-1)-p}{(n-p-1)}$$

Rappel : 

$$R^{2}=\frac{SCE}{SCT}$$

- $CP$ de Mallows : mesure le **biais** introduit par un modèle imbriqué sous-spécifié. Un modèle imbriqué satisfaisant devrai avoir un $Cp$ égal au **nombre de ses

**paramètres**. 

$$Cp=\frac{SCR_{r}}{\sigma^{2}}-(n-2r)$$

Remarques : 

- le *Cp* de Mallows ne peut pas être utilisé sur le modèle complet (def. $Cp=p$)
- dans le cadre d'un modèle **Gaussien**, le *Cp* est équivalent au critères d'Akaike (*AIC*)

- Predicted Residual Sum Squares (*PRESS*) : pertinent dans une optique **prédictive** (LOOCV : Leave One Out Cross-Validation).

$$PRESS=\Sigma^{n}_{i=1}(y-\hat{y}_{(-i)})^{2}$$

Remarque : dans le cadre d'une régression linéaire, on peut simplifier ce calcul :

$$LR \ \Rightarrow \ PRESS=\Sigma^{n}_{i=1}(y-\hat{y}_{(-i)})^{2}=\frac{\Sigma^{n}_{i=1}(y-\hat{y}^{2}}{(1-h_{ii})^{2}}$$

## Indicateurs de stabilité et de multicolinéarité : 

Rappel : 

$$V[B]=\sigma^{2}[X^{T}X]^{-1}$$

Une autocorrélation ayant pour conséquence de réduire certaines qui, faibles dans la matrice d'origine seront trop grandes dans la matrice inverse. On retrouvera cet effet dans une variabilité très forte des paramètres et une mauvaise stabilité du modèle. 

#### Tolérance : 

$$Tol(x_{p})=1-R^{2}_{p}$$

Avec $R^{2}_p$ = coefficient de corrélation pour la régression d'une variable $x_p$ sur les autres variables prédictives.

Propriétés : 

- $Tol\in[0,1]$
- si $Tol=0$, la variable $x_{0}$ est combinaison linéaire des autres

#### Variance Inflation Index (VIF) :

$$VIF(x_{p})=\frac{1}{Tol(x_{p})}=\frac{1}{1-R^{2}_{p}}$$

Propriétés : 

- $VIF \in [1;+\infty]$
- $Tol=1 \ \Rightarrow \ VIF=1$
- on peut réécrire la variance des paramètres par : 

$$\hat\sigma^{2}_{B_{p}}=\frac{\hat{\sigma^{2}}}{\Sigma^{n}_{i=1}(x_{ip}-\bar{x}_{p})^{2}} \times VIF(x_{p})$$

Ordres de grandeur : il est généralement souhaitable de respecter :

- $Tol \geq 0.5$
- $VIF_{x_{p}} \leq 2$

```{r}
SCT=t(Epsilon) %*% Epsilon
# SCE=t(Epsilon) %*% M %*% Epsilon
SCR=t(Epsilon) %*% H %*% Epsilon
```

#### ACP de $(X^{T}X)^{-1}$ :

- *Condition Index* (*Condition Number*) : racine du rapport entre la plus grande et la plus petite des composantes'(ou valeur propre)  dans $(X^{T}X)^{-1}$

- *Proportion of Variation* : indique la proportion de variation expliquée par composantes principales. On évite les fortes proportions (valeurs propres avec un *condition index* élevé).

```{r}
# Estimateur de projection sur le vecteur unitaire moyen :
n=dim(X)[1]
P_unitaire_moyen=matrix(nrow = n,ncol = n,1/n)
# Determination du vecteur de moyennes de Y par projection sur le vecteur unitaire moyen :
moyenne_Y=P_unitaire_moyen%*%Y
# Norme du vecteur Y :
mod_Y = t(Y)%*%Y
# Norme du vecteur Y_barre dans l'espace des observations :
mod_Y_barre = t(Y)%*%P_unitaire_moyen%*%Y
# Norme du vecteur de la somme des carres totale :
SCT=t(Y)%*%(diag(112)-P_unitaire_moyen)%*%Y
SCR=t(Y)%*%(diag(112)-H)%*%Y
SCE=mod_SCR=t(Y)%*%(H-P_unitaire_moyen)%*%Y

# Verification de la complementarite des modules de la moyenne et des carres des residus
# pour former Y :
(as.numeric(mod_Y)-as.numeric(mod_SCT+mod_Y_barre))<10**10
```

```{r}
variables=labels(df)[[2]]
modele = lm(maxO3~T9+T12+T15+Ne9+Ne12+Ne15+Vx9+Vx12+Vx15+maxO3v,data=df)
```

## Thérorèmes et définitions :

### Thérorème de Gauss-Markov :      
Soit $V[ \ b \ ]$, la matrice de covariance de $b$ :

$$V[ \ b \ ]=V[ \ Ay \ ]=A^{T} \ V[ \ y \ ] \ A=\sigma^{2}A \ A^{T}$$
$$=\sigma^{2}(X^{T}X)^{-1}X^{T}X(X^{T}X)^{-1}=\sigma^{2}(X^{T}X)^{-1}$$

Parmi tous les estimateurs lin´eaires non biais´es de β, b est l’estimateur de variance minimale.

Démonstration :

$$\tilde{b}= (A+C) y = Ay+C(X \beta + \epsilon)=Ay+CX \beta + C\epsilon$$

Avec $\tilde{b}$ sans biais :

$$E[ \ \tilde{b} \ ]=E[ \ Ay+CX \beta + C\epsilon \ ]=\beta+CX\beta=\beta(I+CX)=\beta \ \Rightarrow \ CX=0$$

Pour la variance : 

$$V[ \ \tilde{b} \ ]=V[ \ (A+C)Y \ ]=(A+C) \ V[ \ Y \ ] \ (A+C)^{T}=\sigma^{2}(A+C) \ (A+C)^{T}$$
$$=\sigma^{2}(AA^{T}+AC^{T}+CA^{T}+CC^{T})$$

On sait que :       

$$CA^{T}=C(X(X^{T}X)^{-1})=CX(X^{T}X)^{-1} \ \Leftrightarrow \ CX=0 \ \Rightarrow \ CA^{T}=0$$
$$AC^{T}=(CA^{T})^{T}=0$$

Donc : 

$$V[ \ \tilde{b} \ ]=\sigma^{2}AA^{T}+\sigma^{2}CC^{T}$$

Or, $CC^{T}$ étant le carré de $C^{T}$, elle est **symétrique positive**. Donc :

$$V[ \ \tilde{b} \ ]=V[ \ b \ ] \ \ \ \ \Leftrightarrow \ \ \ \ \tilde{b}=b$$

On dit généralement que b est **Best Linear Unbiased Estimator** (BLUE).


### Formes quadratiques (définition) :
Soit $x=\{x_{1},x_{2},...x_{n}\}$ et la matrice $A_{n\times n}$. La forme quadratique de $x$ est le scalaire : 

$$x^{T}Ax=\Sigma_{i=1}^{n}\Sigma_{j=1}^{n}a_{ij}x_{i}x_{j}$$

### Théorème 1 sur les formes quadratiques :
Soit $x=\{x_{1},x_{2},...x_{n}\} \sim \mathcal{N}(0,1)$ et $A$ une matrice de projection orthogonale de $Rang(g)$.

$$Q=x^{T}Ax \sim \chi^{2}_{g}$$

### Théorème de Craig :
Soit $x=\{x_{1},x_{2},...x_{n}\} \sim \mathcal{N}(0,1)$ et $\{A_{1},A_{2}\}$. Alors :

$$\left(Q_{1}=x^{T}A_{1}x\right) \perp \left(Q_{1}=x^{T}A_{1}x\right) \ \ \ SSi \ \ \ A_{1}A_{2}=0$$

### Théorème de Cochran :
Soient $\{Q_{1},...,Q_{k},...,Q_{K} \ |Q_{k}=x^{T}A_{k}x\}$, avec $A_{k}$ une matrice symétrique et $\{\forall i : x \sim \mathcal{N}(0,I_{n})\}$ :

$$\Sigma_{k=1}^{K}Rang(A_{k})=n \ \Rightarrow \ \forall i \neq j : \left(Q_{j} \sim \chi^{2}_{Rang(A_{j})}\right) \perp \left(Q_{j}\sim \chi^{2}_{Rang(A_{j})}\right)$$

**Remarque :** 

$$x \sim \mathcal{N}(0\ ; \ \sigma^{2}I) \ \Rightarrow \ x^{T}x \sim \sigma^{2}\chi^{2}_{n}$$










