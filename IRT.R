# بارگذاری کتابخانه‌ها
library(foreign)    # برای خواندن فایل SPSS
library(ltm)        # مدل‌های IRT ساده مانند 1PL و 2PL
library(mirt)       # مدل‌های IRT پیشرفته‌تر
library(psych)      # تحلیل عاملی و آمار توصیفی
library(ggplot2)    # برای گراف‌ها
library(corrplot)   # برای رسم نقشه‌های حرارتی

# خواندن فایل SPSS 
data <- read.spss(file.choose(), to.data.frame = TRUE, use.value.labels = FALSE)

# استخراج آیتم‌های مربوط به IRT، ستون‌ها A1 تا A10 نام دارند
irt_items <- data[paste0("A", 1:10)]

# تابعی برای اصلاح نمره‌گذاری معکوس (امن نسبت به 0/1 یا 1/2 یا دو مقدار دلخواه)
fix_reversed_binary <- function(x) {
  # تبدیل فاکتور/کاراکتر به عدد (با مهار هشدار)
  if (is.factor(x) || is.character(x)) x_num <- suppressWarnings(as.numeric(as.character(x))) else x_num <- as.numeric(x)
  # اگر قابل تبدیل نبود
  if (all(is.na(x_num) & !is.na(x))) stop("نشد مقدار آیتم را به عدد تبدیل کرد — بررسی کنید A1..A10 چگونه کدگذاری شده‌اند.")
  vals <- sort(unique(na.omit(x_num)))
  if (length(vals) == 0) return(x_num)   # همه NA
  if (length(vals) > 2) stop(paste0("آیتم دارای بیش از دو مقدار متمایز است: ", paste(vals, collapse = ","), ". این کد مناسب آیتم‌های باینری است."))
  # حالت‌های معمول:
  if (all(vals %in% c(0,1))) {
    # مستقیم معکوس کن
    return(ifelse(is.na(x_num), NA, 1 - x_num))
  } else if (all(vals %in% c(1,2))) {
    # ابتدا 1->0, 2->1 و سپس معکوس 
    x01 <- x_num - 1
    return(ifelse(is.na(x_num), NA, 1 - x01))
  } else {
    # حالت عمومی: مقدار کمتر -> 0، بیشتر -> 1، سپس معکوس
    low <- vals[1]; high <- vals[2]
    x01 <- ifelse(x_num == high, 1, ifelse(x_num == low, 0, NA))
    return(ifelse(is.na(x01), NA, 1 - x01))
  }
}

# اعمال اصلاح به همه آیتم‌ها
irt_items_fixed <- as.data.frame(lapply(irt_items, fix_reversed_binary))
names(irt_items_fixed) <- names(irt_items)

# نمایش خلاصه مقادیر قبل و بعد (برای بررسی سریع)
cat("Unique values BEFORE (per item):\n")
print(lapply(irt_items, function(x) sort(unique(na.omit(suppressWarnings(as.numeric(as.character(x))))))))
cat("\nUnique values AFTER (per item) - should be 0 and 1:\n")
print(lapply(irt_items_fixed, function(x) sort(unique(na.omit(x)))))

# توصیف آماری اولیه آیتم‌ها (با داده‌های اصلاح‌شده)
describe(irt_items_fixed)

# محاسبه ماتریس همبستگی تتراکوریک (برای آیتم‌های دو گزینه‌ای)
tetra <- tetrachoric(irt_items_fixed)
cor_matrix <- tetra$rho

# تعداد نمونه‌ها
N <- nrow(data)

# تحلیل عاملی اکتشافی با یک عامل
efa_1factor <- fa(cor_matrix, nfactors = 1, fm = "ml", n.obs = N)
print(efa_1factor$loadings)

# تحلیل عاملی اکتشافی با دو عامل
efa_2factor <- fa(cor_matrix, nfactors = 2, fm = "ml", n.obs = N)
print(efa_2factor$loadings)

# مدل 2PL با پکیج mirt (با داده‌های اصلاح‌شده)
model_mirt <- mirt(irt_items_fixed, 1, itemtype = "2PL")

# بررسی استقلال موضعی با باقیمانده‌های Q3
resid_q3 <- residuals(model_mirt, type = "Q3")

# ذخیره نمودار Q3 به‌صورت Heatmap
png("Q3_heatmap.png", width = 800, height = 600)
corrplot(resid_q3, is.corr = FALSE, method = "color", tl.cex = 0.8)
dev.off()

# باز کردن فایل تصویر
if (Sys.info()["sysname"] == "Windows") {
  try(shell.exec("Q3_heatmap.png"), silent = TRUE)
}

# برازش مدل Rasch (مدل 1PL) با پکیج ltm
rasch_model <- rasch(irt_items_fixed)
summary(rasch_model)

# رسم منحنی ویژگی آیتم‌ها (ICC) برای مدل Rasch
plot(rasch_model, type = "ICC", legend = TRUE)

# برازش مدل Rasch با پکیج mirt (معادل 1PL)
rasch_mirt <- mirt(irt_items_fixed, 1, itemtype = "Rasch")

# رسم تابع آگاهی آزمون برای مدل Rasch 
if (Sys.info()["sysname"] == "Windows") windows()
plot(rasch_mirt, type = "info", main = "Test Information Function - Rasch Model")

# برازش مدل دو پارامتری (2PL) با mirt
model_2pl <- mirt(irt_items_fixed, 1, itemtype = "2PL")

# رسم تابع آگاهی آزمون برای مدل 2PL
if (Sys.info()["sysname"] == "Windows") windows()
plot(model_2pl, type = "info", main = "Test Information Function - 2PL Model")

# مدل 2PL با پکیج ltm (جایگزین دیگر)
model_2PL <- ltm(irt_items_fixed ~ z1)
summary(model_2PL)

# مقایسه مدل Rasch و مدل 2PL از نظر برازش (برای ltm)
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

# باز کردن فایل تصویر 
if (Sys.info()["sysname"] == "Windows") {
  try(shell.exec("ICC_model2PL.png"), silent = TRUE)
}
