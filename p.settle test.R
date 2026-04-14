h.block.mean <- qlogis(0.6)
h.sire <- 0.1
h.dam <- 0.2
h.int <- 0.15

p.hatch <- plogis(h.block.mean + h.sire + h.dam + h.int)
p.hatch


s.gh.block.mean <- qlogis(0.8)
s.gh.sire <- 1
s.gh.dam <- 0.5
s.gh.int <- 1.2

p.settle.gh <- plogis(s.gh.block.mean + s.gh.sire + s.gh.dam + s.gh.int)
p.settle.gh


p.settle <- p.hatch * p.settle.gh
p.settle


p.hatch.wo.sire <- plogis(h.block.mean + h.dam + h.int)
p.hatch.sire <- plogis(h.sire)

p.settle.gh.wo.sire <- plogis(s.gh.block.mean + s.gh.dam + s.gh.int)
p.settle.gh.sire <- plogis(s.gh.sire)

p.settle - (p.hatch.wo.sire * p.settle.gh.wo.sire)
