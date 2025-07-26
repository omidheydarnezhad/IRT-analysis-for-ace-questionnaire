# بارگذاری کتابخانه‌ها
library(foreign)    # برای خواندن فایل SPSS
library(ltm)        # مدل‌های IRT ساده مانند 1PL و 2PL
library(mirt)       # مدل‌های IRT پیشرفته‌تر
library(psych)      # تحلیل عاملی و آمار توصیفی
library(ggplot2)    # برای گراف‌ها
library(corrplot)   # برای رسم نقشه‌های حرارتی

# خواندن فایل SPSS (فایل را از سیستم انتخاب کنید)
data <- read.spss(file.choose(), to.data.frame = TRUE, use.value.labels = FALSE)

# استخراج آیتم‌های مربوط به IRT، فرض بر این است که ستون‌ها A1 تا A10 نام دارند
irt_items <- data[paste0("A", 1:10)]

# توصیف آماری اولیه آیتم‌ها
describe(irt_items)

# محاسبه ماتریس همبستگی تتراکوریک (برای آیتم‌های دو گزینه‌ای)
tetra <- tetrachoric(irt_items)
cor_matrix <- tetra$rho

# تعداد نمونه‌ها
N <- nrow(data)

# تحلیل عاملی اکتشافی با یک عامل
efa_1factor <- fa(cor_matrix, nfactors = 1, fm = "ml", n.obs = N)
print(efa_1factor$loadings)

# تحلیل عاملی اکتشافی با دو عامل
efa_2factor <- fa(cor_matrix, nfactors = 2, fm = "ml", n.obs = N)
print(efa_2factor$loadings)

# مدل 2PL با پکیج mirt
model_mirt <- mirt(irt_items, 1, itemtype = "2PL")

# بررسی استقلال موضعی با باقیمانده‌های Q3
resid_q3 <- residuals(model_mirt, type = "Q3")

# ذخیره نمودار Q3 به‌صورت Heatmap
png("Q3_heatmap.png", width = 800, height = 600)
corrplot(resid_q3, is.corr = FALSE, method = "color", tl.cex = 0.8)
dev.off()

# باز کردن فایل تصویر در ویندوز (اگر در محیط ویندوز هستید)
shell.exec("Q3_heatmap.png")

# برازش مدل Rasch (مدل 1PL) با پکیج ltm
rasch_model <- rasch(irt_items)
summary(rasch_model)

# رسم منحنی ویژگی آیتم‌ها (ICC) برای مدل Rasch
plot(rasch_model, type = "ICC", legend = TRUE)

# برازش مدل Rasch با پکیج mirt (معادل 1PL)
rasch_mirt <- mirt(irt_items, 1, itemtype = "Rasch")

# رسم تابع آگاهی آزمون برای مدل Rasch
windows()  # فقط در ویندوز
plot(rasch_mirt, type = "info", main = "Test Information Function - Rasch Model")

# برازش مدل دو پارامتری (2PL)
model_2pl <- mirt(irt_items, 1, itemtype = "2PL")

# رسم تابع آگاهی آزمون برای مدل 2PL
windows()
plot(model_2pl, type = "info", main = "Test Information Function - 2PL Model")


# مدل 2PL با پکیج ltm (جایگزین دیگر)
model_2PL <- ltm(irt_items ~ z1)
summary(model_2PL)

# مقایسه مدل Rasch و مدل 2PL از نظر برازش
anova(rasch_model, model_2PL)

# رسم منحنی ویژگی آیتم‌ها (ICC) برای مدل 2PL
plot(model_2PL, type = "ICC", legend = TRUE)

# رسم منحنی آگاهی آیتم‌ها (IIC)
plot(model_2PL, type = "IIC", legend = TRUE)

# استخراج ضرایب مدل 2PL
coef(model_2PL)

# ذخیره نمودار ICC برای مدل 2PL در فایل تصویری
png("ICC_model2PL.png", width = 800, height = 600)
plot(model_2PL, type = "ICC", legend = TRUE)
dev.off()

# باز کردن فایل تصویر (اختیاری - فقط در ویندوز)
shell.exec("ICC_model2PL.png")