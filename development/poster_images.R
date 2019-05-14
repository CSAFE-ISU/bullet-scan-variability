library(tidyverse)
library(gridExtra)
library(lme4)
library(bulletxtrctr)


barrel6 <- readRDS("../data/user_aligned_sigs.rda")
b6_long <- barrel6 %>% unnest(sigs_aligned)

sig_by_op <- b6_long %>% 
  filter(Operator != "Marco", Operator != "Jozef") %>%
  filter(unique_id != "BL 6-4", unique_id != "BL 6-3", unique_id != "BL 6-6") %>%
  #filter( unique_id == "BL 6-5") %>% 
  ggplot() + 
  geom_line(aes(x = x_aligned, y = sig2, group = scan_id, color = factor(operator_code)), alpha = 0.7) + 
  #geom_line(aes(x = x_aligned, y = sig1), alpha = 0.8) + 
  theme_bw() + 
  facet_wrap(~unique_id, nrow = 1) + 
  labs(title = "Hamby set Barrel 6", x = "Relative Location",  y = "Signature Height") + 
  scale_color_discrete(name = "Operator")

sig_by_bullet <- b6_long %>% 
  filter(Operator != "Marco", Operator != "Jozef") %>%
  filter(unique_id != "BL 6-4", unique_id != "BL 6-3", unique_id != "BL 6-6") %>%
  #filter( unique_id == "BL 6-5") %>% 
  ggplot() + 
  geom_line(aes(x = x_aligned, y = sig2, group = scan_id, color = factor(bullet_code)), alpha = 0.7) + 
  #geom_line(aes(x = x_aligned, y = sig1), alpha = 0.8) + 
  theme_bw() + 
  facet_wrap(~unique_id, nrow = 1) + 
  labs(title = "Hamby set Barrel 6", x = "Relative Location",  y = "Signature Height") + 
  scale_color_discrete(name = "Bullet")





pairwise_data <- readRDS("data/user_meas_pairwise.rda")
bullet_gt <- read_csv("data/user_bullets_gt.csv")


bullet_s1 <- bullet_gt %>% 
  mutate(Operator_S1 = Operator, 
         Barrel_S1 = Barrel, 
         Bullet_S1 = Bullet, 
         Round_S1 = Round, 
         bullet_id_S1 = bullet_id, 
         bullet_code_S1 = bullet_code) %>%
  select(Operator_S1, Barrel_S1, Bullet_S1, Round_S1, bullet_id_S1, bullet_code_S1) 

bullet_s2 <- bullet_gt %>% 
  mutate(Operator_S2 = Operator, 
         Barrel_S2 = Barrel, 
         Bullet_S2 = Bullet, 
         Round_S2 = Round, 
         bullet_id_S2 = bullet_id, 
         bullet_code_S2 = bullet_code) %>%
  select(Operator_S2, Barrel_S2, Bullet_S2, Round_S2, bullet_id_S2, bullet_code_S2) 

pairwise_data <- left_join(pairwise_data, bullet_s1)
pairwise_data <- left_join(pairwise_data, bullet_s2)




pairwise_data <- pairwise_data %>% 
  filter(unique_id != "BL 6-4", unique_id != "BL 6-3", unique_id != "BL 6-6") %>%
  mutate(operator_truth = ifelse(Operator_S1 == Operator_S2, 
                                 "Same Operator", "Diff Operator"), 
         bullet_truth = ifelse(bullet_code_S1 == bullet_code_S2,
                               "Same Bullet", "Diff Bullet"),
         machine_truth = ifelse(Machine_S1 == Machine_S2,
                                "Same Machine", "Diff Machine"),
         round_truth = ifelse(Round_S1 == Round_S2,
                              "Same Round", "Diff Round"))
head(pairwise_data)



pairwise_data %>%
  filter(ground_truth == "Same Source") %>%
  filter(Operator_S1 != "Marco", Operator_S1 != "Jozef") %>%
  filter(Operator_S2 != "Marco", Operator_S2 != "Jozef") %>%
  filter(Barrel_S1 == "Barrel 6") %>%
  filter(unique_id_S1 != "BL 6-4", unique_id_S1 != "BL 6-3", unique_id_S1 != "BL 6-6") %>%
  #filter(rfscore > 0.55) %>%
  ggplot() + 
  geom_density(aes(x = rfscore, fill = factor(bullet_truth)), alpha = 0.6) + 
  theme_bw() + 
  facet_wrap(~unique_id_S1, nrow = 3) + 
  labs(x = "Random Forest Score", y = "Density") + 
  scale_fill_discrete(name = "Bullet Pairing")