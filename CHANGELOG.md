# Changelog

## 1.0.0

- The four stand-alone repositories RDE26-SOP, RDE26-CSOP, RDE26-MOP and RDE26-CMOP are combined in this repository, one directory per track (`RDE26_SOP`, `RDE26_CSOP`, `RDE26_MOP`, `RDE26_CMOP`), in the layout of the RDEx series. The algorithm code is unchanged from the heads of the four repositories on 2026-07-09 (RDE26-SOP 23732ed, RDE26-CSOP 46e5698, RDE26-MOP 1fdedd2, RDE26-CMOP dc2c1f6); only the folder locations changed.
- Each track builds and runs from its new location (checked with short runs of every track).
- Added: top-level and per-track READMEs in English, Chinese and Japanese, `LICENSE` (MIT), `CITATION.cff`, `THIRD_PARTY_NOTICES.md`, and this file.
- The READMEs of the four stand-alone repositories (navigation tables, repeated paper lists) are replaced by the READMEs above. Scores are now given once, in the top-level README; the BC-SOP score follows version V6 of the organizers' slides (the earlier READMEs gave the value of an earlier version of the slides).
- `.gitignore` of the four repositories merged, with the paths adapted to the new layout.
