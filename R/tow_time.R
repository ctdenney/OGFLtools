# Tow time helpers
#
# Field datasheet times appear as "HHMM" (e.g. "0915", "915") or "HH:MM", with
# stray separators/whitespace. parse_tow_time() returns minutes since midnight
# (for ordering/comparison); format_tow_time() returns a canonical "HH:MM"
# string. Both yield NA for blank, malformed, or out-of-range (> 23:59) values.

parse_tow_time <- function(x) {
  d  <- gsub('[^0-9]', '', as.character(x))
  n  <- nchar(d)
  hh <- ifelse(n == 4, substr(d, 1, 2),
        ifelse(n == 3, substr(d, 1, 1), NA_character_))
  mm <- ifelse(n %in% c(3, 4), substr(d, n - 1, n), NA_character_)
  hh <- suppressWarnings(as.integer(hh))
  mm <- suppressWarnings(as.integer(mm))
  out <- hh * 60L + mm
  out[is.na(hh) | is.na(mm) | hh > 23 | mm > 59] <- NA_integer_
  out
}

format_tow_time <- function(x) {
  m   <- parse_tow_time(x)
  out <- rep(NA_character_, length(m))
  ok  <- !is.na(m)
  out[ok] <- sprintf('%02d:%02d', m[ok] %/% 60L, m[ok] %% 60L)
  out
}
