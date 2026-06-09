### Version 0.0.1 (June 2025)
- Initial Release

### Version 0.0.2 (May 2026)
- Data can be simulated in parallel for the quantile model.
- Bootstrap methods are added for reliable inference.
- The Encompassing test now uses bootstrap because compared models can be misspecified
- Optimal instruments are added (see the `qpeer.optimal.inst` function) but requires iid errors (see Herstad, Houndetoungan, and Shin, 2026).
- Some notation changes have been made to ensure consistency with the new version of the paper (see Herstad, Houndetoungan, and Shin, 2026). In particular, the list of network matrices is now denoted by `Glist` rather than `A`. 