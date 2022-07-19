# age_con = function(ages,error_type,error_lim,read_col, age_col, readid_col,fishid_col,reader_col) {
#   if(error_type == 'cv') {
#     a = with(ages, aggregate(age, list(reader_col,fishid_col,readid_col), FUN = max))
#     a_by = by(
#       a,
#       INDICES = a[,fishid_col],
#       FUN = funtion(x) {
#         x$R = length(x[,readid_col])
#         return(x)
#       }
#     )
#     a_comb = do.call(rbind, a_by)
#     return(a)
#   }
# }
#
# test = age_con(
#   ages = pass_reads, error_type = 'ape',error_lim = 10,
#   read_col = 'obs',age_col = 'age',readid_col = 'readid',fishid_col = 'fishid'
# )
#
# test = with(pass_reads, by(fishid, readid, obs, FUN = max))
#
# test = with(pass_reads, aggregate(age, list(obs,fishid, readid), FUN = max))
#
# test_by = by(
#   pass_reads,
#   INDICES = pass_reads[,c('obs','fishid','readid')],
#   FUN = function(x) {
#     data.frame(
#       age = max(age)
#     )
#   }
# )
# test =
#
#   names(test) = c('obs','fishid','readid','age')
#
# test_by = by(
#   test,
#   INDICES = test[,'fishid'],
#   FUN = function(x) {
#     x$R = length(x[,'readid'])
#     return(x)
#   }
# )
#
# test_comb = do.call(rbind, test_by)
#
# error = with(test_comb, aggregate())
