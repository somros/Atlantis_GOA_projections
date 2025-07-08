# build hcr shapes for proposal plot

# make a data frame
df <- data.frame("type" = c(rep("SQ", 1001),
                            rep("SQ_SSL", 1001),
                            rep("SQ_B50", 1001),
                            rep("SQ_SSL_B50", 1001)),
                 "ssl" = rep(c(rep(0,1001), rep(1,1001)), 2),
                 "ssb" = rep(seq(0,1,0.001), 4),
                 "alpha" = 0.05,
                 "btarg" = c(rep(0.4,2002),
                             rep(0.5,2002)))


make_hcr <- function(ssb,ssl,alpha,btarg){
  
  if(ssb/btarg > 1) 
  {f <- 1}
  else if(ssb/btarg > alpha & ssb/btarg <= 1)
  { f <- 1 * (ssb/btarg - alpha) / (1 - alpha)}
  else 
  { f <- 0}
  
  if(ssl==1){
    if(ssb/btarg < 0.5)
    {f <- 0}
  }
  
  return(f)
  
}

df <- df %>%
  rowwise() %>%
  mutate(f = make_hcr(ssb, ssl, alpha, btarg)) %>%
  ungroup()

# plot in ggplot
p <- df %>%
  ggplot(aes(x = ssb, y = f))+
  geom_line(linewidth = 1)+
  theme_bw()+
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey")+
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "blue")+
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "red")+
  scale_x_continuous(breaks = seq(0,1,0.1))+
  labs(x = "B/B0", y = "F/Ftarg")+
  facet_wrap(~type)
p

ggsave("hcr_nprb.png",p,width = 8, height = 4)
