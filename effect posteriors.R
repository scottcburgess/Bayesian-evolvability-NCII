pdf('sire_effects.pdf')
plot(post, vars = 'sire.eff')
dev.off()

pdf('dam_effects.pdf')
plot(post, vars = 'dam.eff')
dev.off()

pdf('int_effects.pdf')
plot(post, vars = 'int.eff')
dev.off()

pdf('means.pdf')
plot(post, vars = c('block.mean', 'sire.mean', 'dam.mean', 'int.mean'))
dev.off()

