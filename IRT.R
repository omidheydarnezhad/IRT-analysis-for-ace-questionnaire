# -----------------------------
#  بارگذاری کتابخانه‌ها
# -----------------------------
library(foreign)
library(ltm)
library(mirt)
library(psych)
library(ggplot2)
library(corrplot)

# -----------------------------
#  خواندن فایل SPSS
# -----------------------------
data <- read.spss(file.choose(), to.data.frame = TRUE, use.value.labels = FALSE)

# -----------------------------
#  استخراج آیتم‌های IRT 
#  (نام آیتم‌ها اکنون a1 تا a10 هستند)
# -----------------------------
irt_items <- data[paste0("a", 1:10)]

# تبدیل فاکتور/کاراکتر به عدد برای جلوگیری از خطا
irt_items <- as.data.frame(lapply(irt_items, function(x) {
  if (is.factor(x) || is.character(x)) as.numeric(as.character(x)) else x
}))

# -----------------------------
#  توصیف آماری آیتم‌ها
# -----------------------------
describe(irt_items)

# -----------------------------
#  ماتریس همبستگی تتراکوریک
# -----------------------------
tetra <- tetrachoric(irt_items)
cor_matrix <- tetra$rho

# تحلیل موازی برای EFA
fa.parallel(tetra$rho,
            n.obs = nrow(irt_items),
            fa = "fa",
            fm = "minres")
# -----------------------------
#  تحلیل عاملی یک‌عاملی
# -----------------------------
N <- nrow(data)
efa_1factor <- fa(cor_matrix, nfactors = 1, fm = "ml", n.obs = N)
print(efa_1factor$loadings)

# -----------------------------
#  تحلیل عاملی دو‌عاملی
# -----------------------------
efa_2factor <- fa(cor_matrix, nfactors = 2, fm = "ml", n.obs = N)
print(efa_2factor$loadings)

# -----------------------------
#  مدل 2PL با mirt
# -----------------------------
model_mirt <- mirt(irt_items, 1, itemtype = "2PL")

# -----------------------------
#  Q3 Residuals (Local Dependence)
# -----------------------------
resid_q3 <- residuals(model_mirt, type = "Q3")

png("Q3_heatmap.png", width = 800, height = 600)
corrplot(resid_q3, is.corr = FALSE, method = "color", tl.cex = 0.8)
dev.off()

if (Sys.info()["sysname"] == "Windows") shell.exec("Q3_heatmap.png")

# -----------------------------
#  مدل Rasch (1PL) با ltm
# -----------------------------
rasch_model <- rasch(irt_items)
summary(rasch_model)

plot(rasch_model, type = "ICC", legend = TRUE)

# -----------------------------
#  مدل Rasch با mirt
# -----------------------------
rasch_mirt <- mirt(irt_items, 1, itemtype = "Rasch")

if (Sys.info()["sysname"] == "Windows") windows()
plot(rasch_mirt, type = "info", main = "Test Information Function - Rasch Model")

# -----------------------------
#  مدل 2PL با mirt
# -----------------------------
model_2pl <- mirt(irt_items, 1, itemtype = "2PL")

if (Sys.info()["sysname"] == "Windows") windows()
plot(model_2pl, type = "info", main = "Test Information Function - 2PL Model")

# -----------------------------
#  مدل 2PL با ltm
# -----------------------------
model_2PL <- ltm(irt_items ~ z1)
summary(model_2PL)

anova(rasch_model, model_2PL)

plot(model_2PL, type = "ICC", legend = TRUE)
plot(model_2PL, type = "IIC", legend = TRUE)

# -----------------------------
#  ضرایب مدل 2PL
# -----------------------------
coef(model_2PL)

# -----------------------------
#  ذخیره ICC مدل 2PL
# -----------------------------
png("ICC_model2PL.png", width = 800, height = 600)
plot(model_2PL, type = "ICC", legend = TRUE)
dev.off()

if (Sys.info()["sysname"] == "Windows") shell.exec("ICC_model2PL.png")
# محاسبه ماتریس همبستگی تتراکوریک
tetra <- tetrachoric(irt_items)
cor_matrix <- tetra$rho

# آزمون KMO
KMO(cor_matrix)

# آزمون Bartlett
cortest.bartlett(cor_matrix, n = nrow(irt_items))

# مدل‌ها فرضی: rasch_mirt (1PL) و model_2pl (2PL) از قبل برازش شده‌اند
anova(rasch_mirt, model_2pl)

model_2pl <- mirt(irt_items, 1, itemtype = "2PL")
