# 准备 ---------------------------------------------------------------

## library
library(GEOquery)
library(SummarizedExperiment)
library(QCnormSE)
library(dplyr)

rm(list = ls())
## dir_create
fs::dir_create("./2-GEO下载")

## 全局参数
gse_id <- "60670"

# 01-预处理-------------------------------------------------

## GEO数据下载
accession_string <- paste0("GSE", gse_id)
se <- QCnormSE::get_geo_data(
  accession = accession_string,
  # platform = "GPL97", # 多个平台再填写， #nolint: commented_code_linter.
  temp.dir = "./2-GEO下载"
)
se
save(se, file = paste0("./2-GEO下载/", gse_id, ".rda"))
load(paste0("./2-GEO下载/", gse_id, ".rda"))
## 肉眼看摘要
row_dat <- rowData(se) %>% as.data.frame()
col_dat <- colData(se) %>% as.data.frame()

# 02-清洗基因名-------------------------------------------------

## 发现是GB_ACC：注释为基因名
if ("GB_ACC" %in% colnames(row_dat)) {
  library(org.Hs.eg.db)
  columns(org.Hs.eg.db)
  a <- AnnotationDbi::select(org.Hs.eg.db,
    keys = row_dat$GB_ACC, keytype = "ACCNUM",
    columns = c("SYMBOL")
  )
  rowData(se)$gene_id <- a$SYMBOL
  se <- se[!is.na((rowData(se)$gene_id)), ]
}

## 发现Gene.symbol：清洗后添加gene_id
if ("Gene.symbol" %in% colnames(row_dat)) {
  keep <-
    row_dat %>%
    dplyr::filter(!grepl("///", .$Gene.symbol)) %>%
    ## 去除基因名中有///的
    dplyr::filter(.$Gene.symbol != "") %>%
    ## 去除基因名为空的
    dplyr::filter(!grepl("---", .$Gene.symbol)) %>%
    ## 去除基因名中有---的
    dplyr::rename(gene_id = "Gene.symbol")

  se <- se[rownames(se) %in% keep$ID, ]
  rowData(se)$gene_id <- keep$gene_id
}

if (FALSE) {
  keep <-
    row_dat %>%
    # dplyr::mutate(gene_id = stringr::str_extract(.$gene_assignment, "(?<= // )[^/]+")) %>% # nolint
    dplyr::mutate(gene_id = stringr::str_split(.$gene_assignment, " // ", simplify = TRUE)[, 2]) %>% # nolint
    dplyr::filter(!grepl("---", .$gene_id))

  se <- se[rownames(se) %in% keep$ID, ]
  rowData(se)$gene_id <- keep$gene_id
}

### 整理分组
se$title

# ## 过滤不需要的样本
# se <- se[, !grepl("0h", col_dat$title)]
# se

se$group <- case_when(
  grepl("Sham", se$title) ~ "Normal",
  # grepl("GP", se$title) ~ "Disorder",
  TRUE ~ "Disorder"
)
table(se$group)

se <- se[, se$group %in% c("Normal", "Disorder")]

## 注释
sum(duplicated(rowData(se)$gene_id)) # 重复基因计数
se <- se[!duplicated(rowData(se)$gene_id), ]

rownames(se) <- rowData(se)$gene_id
se

# 保存-------------------------------------------------

## 提取数据
exprs <- assays(se)$exprs
group <- colData(se) %>%
  as.data.frame() %>%
  dplyr::select("geo_accession", "group")
colnames(group) <- c("sample", "group")

## 查看
range(exprs)
sum(is.na(exprs))

## 文件路径
exprs_file_path <- sprintf("./2-GEO下载/exprs_%s.csv", gse_id)
group_file_path <- sprintf("./2-GEO下载/group_%s.csv", gse_id)

## 使用构建的文件路径来保存CSV文件
write.csv(exprs, exprs_file_path, quote = FALSE)

write.csv(group, group_file_path, row.names = FALSE, quote = FALSE)
