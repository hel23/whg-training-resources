---
sidebar_position: 16
---

# Appendix: visualisation

[Up to table of contents](README.md)

It would be very nice to plot some of the values we've computed! For example - gene length versus
number of exons, or gene length versus number of transcripts, or genome size versus number of
genes, and so on.

Data visualisation is really beyond the scope of this tutorial. However, as a quick example, here is how you might plot some of
these values in R using [ggplot2](https://ggplot2.tidyverse.org). (You need to have created the
[SQL summaries](Counting_genes_2.md#a-sqlite-approach) first.)

```R
library( RSQLite )
library( ggplot2 )
db = DBI::dbConnect( RSQLite::SQLite(), "genes.sqlite" )
data = dbGetQuery( db, "SELECT * FROM gene_summary_view" )
data$length = data$end - data$start + 1

p = (
   ggplot( data = data )
   + geom_point( aes( x = end - start + 1, y = average_number_of_exons, col = analysis ))
   + xlab( "Gene length" )
   + ylab( "Number of exons" )
   + facet_wrap( ~analysis )
)

ggsave( p, file = "gene_length_versus_number_of_exons.pdf", width = 8, height = 5 )
```

![img](images/visualisation_example.png)

In python, you could instead try using  [matplotlib](https://matplotlib.org). I'd like to
include a matplotlib example here but don't have one - let me know if you write one!

## Conclusion

That's the end of this tutorial!  You can read some [closing thoughts](Closing_thoughts.md).
