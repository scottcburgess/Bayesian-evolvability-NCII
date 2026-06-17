pdf('Model_outputs/sire_effects.pdf')
plot(post, vars = 'sire.eff')
dev.off()

pdf('Model_outputs/dam_effects.pdf')
plot(post, vars = 'dam.eff')
dev.off()

pdf('Model_outputs/int_effects.pdf')
plot(post, vars = 'int.eff')
dev.off()

pdf('Model_outputs/means.pdf')
plot(post, vars = c('block.mean', 'sire.mean', 'dam.mean', 'int.mean'))
dev.off()

