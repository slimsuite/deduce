from typing import Any, Dict
from deduce_uces.mappers.bowtie_mapper import BowtieMapper
from deduce_uces.mappers.minimap_mapper import MinimapMapper

MAPPING_STRATEGIES: Dict[str, Any] = {"minimap": MinimapMapper, "bowtie": BowtieMapper}
