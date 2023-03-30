\dontshow{
if (requireNamespace("CCA", quietly = TRUE)) \{
}
data(nutrimouse, package = "CCA")

x <- nutrimouse$gene[ , 1:5]
y <- nutrimouse$lipid
cc <- cancor(x, y)

# Re-compute explained correlation
ac <- acor(x, cc$xcoef, y, cc$ycoef)

# Results should agree
print(cc$cor)
print(ac$cor)
\dontshow{\}}
