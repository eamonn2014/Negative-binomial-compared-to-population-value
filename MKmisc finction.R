

require(MKpower)
# ssize.ratio = n1/n0


# Approach 1 underestimates the sample size because it ignores overdispersion.
# 
# Approach 2 provides a more accurate estimate by accounting for overdispersion.
# 
# Approach 3 is the most conservative, yielding a larger sample size to account for additional uncertainty in the dispersion parameter.


# The function should calculate the required sample size for the treatment group (n1) 
# based on the provided parameters. Since ssize.ratio = 1e6, the reference group
# is effectively treated as a known population, and the result will approximate a one-sample test.

# 
# n = control sample size
# ssize.ratio = n1/n, set n to massive =1/1e6
# mu0 = control mean =1 
# interest in n1 the ssize for grp1 
# theya = 1/k

# 59
power.nb.test(n = NULL, mu0=1, mu1=.58, duration = 1, theta=1/1.78, ssize.ratio = 1/1e7,  ##n1/n
                      sig.level = 0.01, power = .8, alternative = c("one.sided"),
                       approach = 1)

power.nb.test(n = NULL, mu0=1, RR=.58, duration = 1, theta=1/1.78, ssize.ratio = 1/1e6,  ##n1/n
              sig.level = 0.01, power = .8, alternative = c("one.sided"),
              approach = 1)


# 78
power.nb.test(n = NULL, mu0=1, mu1=.58, duration = 1, theta=1/1.78, ssize.ratio = 1/1e6,
              sig.level = 0.01, power = .8, alternative = c("one.sided"),
              approach = 2)


# 59
power.nb.test(n = NULL, mu0=1, mu1=.58, duration = 1, theta=1/1.78, ssize.ratio = 1/1e6,
              sig.level = 0.01, power = .8, alternative = c("one.sided"),
              approach = 3)



# these funbctions dont agree 
require(skewsamp)

n_negbinom(
  mean0 =71.4,
  effect = .3,
  dispersion0 =.33,
  dispersion1 = .33,
  alpha = 0.05,
  power = 0.9,
  q = 0.5,
  link = c("identity"),
  two_sided = FALSE
)

# this matches the paper!!!! skewed distributins paper example, I enter k for theta??
power.nb.test(n = NULL, mu0=71.4, mu1=50, duration = 1, theta=.33, ssize.ratio = 1,
              sig.level = 0.05, power = .9, alternative = c("two.sided"),
              approach = 3)












# these funbctions dont agree 
require(skewsamp)

n_negbinom(
  mean0 =1,
  effect =1-.58/1,
  dispersion0 =1/1.78,
  dispersion1 = 1/1.78,
  alpha = 0.01,
  power = 0.8,
  q = 0.5,
  link = c("identity"),
  two_sided = FALSE
)



power.nb.test(n = NULL, mu0=1, mu1=.58, duration = 1, theta=1/1.78, ssize.ratio = 1,
              sig.level = 0.01, power = .8, alternative = c("one.sided"),
              approach = 3)
