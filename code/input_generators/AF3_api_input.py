# custom class for the creation of input files for the AF3 API

from __future__ import annotations

import copy
import json
import re
import warnings
from pathlib import Path
from typing import Any, Mapping


class _FieldSpec:
	def __init__(self, value_types: type | tuple[type, ...], item_type: type | tuple[type, ...] | None = None, nested_class: type | None = None) -> None:
		if not isinstance(value_types, tuple):
			value_types = (value_types,)
		self.value_types = value_types
		self.item_type = item_type
		self.nested_class = nested_class

	def validate(self, value: Any) -> bool:
		if value is None:
			return True

		if not isinstance(value, self.value_types):
			return False

		if self.nested_class is not None:
			for item in value:
				if not isinstance(item, self.nested_class) and not isinstance(item, Mapping):
					return False

		return self._validate_items(value)

	def coerce(self, value: Any) -> Any:
		if value is None:
			return None

		if self.nested_class is not None and isinstance(value, list):
			coerced_items = []
			for item in value:
				if isinstance(item, self.nested_class):
					coerced_items.append(item)
				elif isinstance(item, Mapping):
					coerced_items.append(self.nested_class(item))
				else:
					coerced_items.append(item)
			return coerced_items

		return value

	def _validate_items(self, value: Any) -> bool:
		if self.item_type is None or not isinstance(value, list):
			return True

		for item in value:
			if self.item_type is int:
				if not isinstance(item, int) or isinstance(item, bool):
					return False
			elif not isinstance(item, self.item_type):
				return False

		return True


class _AF3JSONBase:
	_JSON_KEY: str | None = None
	_FIELD_SPECS: dict[str, _FieldSpec] = {}
	_DEFAULTS: dict[str, Any] = {}

	def __init__(self, data: Mapping[str, Any] | None = None, **kwargs: Any) -> None:
		for field_name, default_value in self._DEFAULTS.items():
			setattr(self, field_name, copy.deepcopy(default_value))

		self.update(data=data, **kwargs)

	def update(self, data: Mapping[str, Any] | None = None, **kwargs: Any) -> _AF3JSONBase:
		payload: dict[str, Any] = {}
		if data is not None:
			payload.update(dict(data))
		payload.update(kwargs)

		for key, value in payload.items():
			self._set_field(key, value)

		self._post_update()
		return self

	@classmethod
	def from_json_file(cls, json_path: str | Path) -> _AF3JSONBase:
		path = Path(json_path)
		with path.open("r", encoding="utf-8") as handle:
			loaded = json.load(handle)

		if not isinstance(loaded, dict):
			warnings.warn(
				f"Expected a JSON object in {path}, got {type(loaded).__name__}. No fields were loaded.",
				stacklevel=2,
			)
			return cls()

		if cls._JSON_KEY and cls._JSON_KEY in loaded and isinstance(loaded[cls._JSON_KEY], dict):
			loaded = loaded[cls._JSON_KEY]

		return cls(loaded)

	def to_payload(self) -> dict[str, Any]:
		data: dict[str, Any] = {}
		for field_name in self._FIELD_SPECS:
			data[field_name] = self._serialize_value(getattr(self, field_name))
		return {key: value for key, value in data.items() if value is not None}

	def to_dict(self) -> dict[str, Any]:
		payload = self.to_payload()
		if self._JSON_KEY is None:
			return payload
		return {self._JSON_KEY: payload}

	def write_json(self, path: str | Path, prefix: str = "") -> Path:
		target = Path(path)
		if target.suffix.lower() != ".json":
			target.mkdir(parents=True, exist_ok=True)
			file_name = f"{prefix}{self._safe_file_stem()}.json"
			target = target / file_name
		else:
			target.parent.mkdir(parents=True, exist_ok=True)

		with target.open("w", encoding="utf-8") as handle:
			json.dump(self.to_dict(), handle, indent=2)

		return target

	def _set_field(self, key: str, value: Any) -> None:
		spec = self._FIELD_SPECS.get(key)
		if spec is None:
			warnings.warn(f"Ignoring unknown AF3 input field: {key!r}", stacklevel=2)
			return

		if not spec.validate(value):
			warnings.warn(
				f"Ignoring AF3 input field {key!r} because the value has an unexpected type: {type(value).__name__}.",
				stacklevel=2,
			)
			return

		setattr(self, key, spec.coerce(value))

	def _post_update(self) -> None:
		return

	def _serialize_value(self, value: Any) -> Any:
		if isinstance(value, _AF3JSONBase):
			return value.to_dict() if value._JSON_KEY is not None else value.to_payload()

		if isinstance(value, list):
			return [self._serialize_value(item) for item in value]

		return value

	def _safe_file_stem(self) -> str:
		for attr_name in ("name", "id"):
			if not hasattr(self, attr_name):
				continue

			value = getattr(self, attr_name)
			if value is None:
				continue

			if isinstance(value, list):
				value = "_".join(str(item) for item in value)

			text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._-")
			if text:
				return text

		return self.__class__.__name__.lower()


class ProteinModification(_AF3JSONBase):
	_JSON_KEY = None
	_FIELD_SPECS = {
		"ptmType": _FieldSpec(str),
		"ptmPosition": _FieldSpec(int),
	}
	_DEFAULTS = {
		"ptmType": None,
		"ptmPosition": None,
	}


class RnaModification(_AF3JSONBase):
	_JSON_KEY = None
	_FIELD_SPECS = {
		"modificationType": _FieldSpec(str),
		"basePosition": _FieldSpec(int),
	}
	_DEFAULTS = {
		"modificationType": None,
		"basePosition": None,
	}


class DnaModification(_AF3JSONBase):
	_JSON_KEY = None
	_FIELD_SPECS = {
		"modificationType": _FieldSpec(str),
		"basePosition": _FieldSpec(int),
	}
	_DEFAULTS = {
		"modificationType": None,
		"basePosition": None,
	}


class Template(_AF3JSONBase):
	_JSON_KEY = None
	_FIELD_SPECS = {
		"mmcif": _FieldSpec(str),
		"mmcifPath": _FieldSpec(str),
		"queryIndices": _FieldSpec(list, item_type=int),
		"templateIndices": _FieldSpec(list, item_type=int),
	}
	_DEFAULTS = {
		"mmcif": None,
		"mmcifPath": None,
		"queryIndices": None,
		"templateIndices": None,
	}

	def _post_update(self) -> None:
		if self.mmcif is not None and self.mmcifPath is not None:
			warnings.warn("Template allows either mmcif or mmcifPath, not both.", stacklevel=2)


class Protein(_AF3JSONBase):
	_JSON_KEY = "protein"
	_FIELD_SPECS = {
		"id": _FieldSpec((str, list), item_type=str),
		"sequence": _FieldSpec(str),
		"modifications": _FieldSpec(list, nested_class=ProteinModification),
		"description": _FieldSpec(str),
		"unpairedMsa": _FieldSpec(str),
		"unpairedMsaPath": _FieldSpec(str),
		"pairedMsa": _FieldSpec(str),
		"pairedMsaPath": _FieldSpec(str),
		"templates": _FieldSpec(list, nested_class=Template),
	}
	_DEFAULTS = {
		"id": None,
		"sequence": None,
		"modifications": [],
		"description": None,
		"unpairedMsa": None,
		"unpairedMsaPath": None,
		"pairedMsa": None,
		"pairedMsaPath": None,
		"templates": [],
	}

	def _post_update(self) -> None:
		if self.unpairedMsa is not None and self.unpairedMsaPath is not None:
			warnings.warn("Protein allows either unpairedMsa or unpairedMsaPath, not both.", stacklevel=2)
		if self.pairedMsa is not None and self.pairedMsaPath is not None:
			warnings.warn("Protein allows either pairedMsa or pairedMsaPath, not both.", stacklevel=2)
		if (self.unpairedMsa is None) != (self.pairedMsa is None):
			warnings.warn("Protein unpairedMsa and pairedMsa must both be set or both be unset.", stacklevel=2)


class RNA(_AF3JSONBase):
	_JSON_KEY = "rna"
	_FIELD_SPECS = {
		"id": _FieldSpec((str, list), item_type=str),
		"sequence": _FieldSpec(str),
		"modifications": _FieldSpec(list, nested_class=RnaModification),
		"description": _FieldSpec(str),
		"unpairedMsa": _FieldSpec(str),
		"unpairedMsaPath": _FieldSpec(str),
	}
	_DEFAULTS = {
		"id": None,
		"sequence": None,
		"modifications": [],
		"description": None,
		"unpairedMsa": None,
		"unpairedMsaPath": None,
	}

	def _post_update(self) -> None:
		if self.unpairedMsa is not None and self.unpairedMsaPath is not None:
			warnings.warn("RNA allows either unpairedMsa or unpairedMsaPath, not both.", stacklevel=2)


class DNA(_AF3JSONBase):
	_JSON_KEY = "dna"
	_FIELD_SPECS = {
		"id": _FieldSpec((str, list), item_type=str),
		"sequence": _FieldSpec(str),
		"modifications": _FieldSpec(list, nested_class=DnaModification),
		"description": _FieldSpec(str),
	}
	_DEFAULTS = {
		"id": None,
		"sequence": None,
		"modifications": [],
		"description": None,
	}


class Ligand(_AF3JSONBase):
	_JSON_KEY = "ligand"
	_FIELD_SPECS = {
		"id": _FieldSpec((str, list), item_type=str),
		"ccdCodes": _FieldSpec(list, item_type=str),
		"smiles": _FieldSpec(str),
		"description": _FieldSpec(str),
	}
	_DEFAULTS = {
		"id": None,
		"ccdCodes": None,
		"smiles": None,
		"description": None,
	}

	def _post_update(self) -> None:
		if self.ccdCodes is not None and self.smiles is not None:
			warnings.warn("Ligand allows either ccdCodes or smiles, not both.", stacklevel=2)


class AF3_input(_AF3JSONBase):
	_JSON_KEY = None
	_FIELD_SPECS = {
		"name": _FieldSpec(str),
		"modelSeeds": _FieldSpec(list, item_type=int),
		"sequences": _FieldSpec(list),
		"bondedAtomPairs": _FieldSpec(list),
		"userCCD": _FieldSpec(str),
		"userCCDPath": _FieldSpec(str),
		"dialect": _FieldSpec(str),
		"version": _FieldSpec(int),
	}
	_DEFAULTS = {
		"name": None,
		"modelSeeds": None,
		"sequences": [],
		"bondedAtomPairs": None,
		"userCCD": None,
		"userCCDPath": None,
		"dialect": "alphafold3",
		"version": 4,
	}

	def set_seed(self, seed: int) -> AF3_input:
		if not isinstance(seed, int) or isinstance(seed, bool):
			raise TypeError("seed must be an int")

		self.modelSeeds = [seed]
		return self

	def set_model_seeds(self, model_seeds: list[int]) -> AF3_input:
		if not isinstance(model_seeds, list):
			raise TypeError("model_seeds must be a list of ints")
		if not model_seeds:
			raise ValueError("model_seeds must contain at least one seed")
		if any((not isinstance(seed, int) or isinstance(seed, bool)) for seed in model_seeds):
			raise TypeError("all model seeds must be ints")

		self.modelSeeds = list(model_seeds)
		return self

	def choose_model_seeds(self, model_seeds: list[int]) -> AF3_input:
		return self.set_model_seeds(model_seeds)

	def add_sequence_from_type(self, sequence_type: str, sequence_object: Protein | RNA | DNA | Ligand) -> AF3_input:
		sequence_type_normalized = sequence_type.lower()
		sequence_class_map = {
			"protein": Protein,
			"rna": RNA,
			"dna": DNA,
			"ligand": Ligand,
		}

		sequence_class = sequence_class_map.get(sequence_type_normalized)
		if sequence_class is None:
			raise ValueError("sequence_type must be one of: protein, rna, dna, ligand")
		if not isinstance(sequence_object, sequence_class):
			raise TypeError(f"sequence_object must be an instance of {sequence_class.__name__}")

		self.sequences.append(sequence_object.to_dict())
		return self

	def add_sequence(self, sequence_object: Protein | RNA | DNA | Ligand) -> AF3_input:
		sequence_class_map = (
			(Protein, "protein"),
			(RNA, "rna"),
			(DNA, "dna"),
			(Ligand, "ligand"),
		)

		for sequence_class, sequence_type in sequence_class_map:
			if isinstance(sequence_object, sequence_class):
				return self.add_sequence_from_type(sequence_type, sequence_object)

		raise TypeError("sequence_object must be an instance of Protein, RNA, DNA, or Ligand")

	def _safe_file_stem(self) -> str:
		if self.name:
			text = re.sub(r"[^A-Za-z0-9._-]+", "_", str(self.name)).strip("._-")
			if text:
				return text
		return "af3_input"


__all__ = [
	"AF3_input",
	"Protein",
	"RNA",
	"DNA",
	"Ligand",
	"ProteinModification",
	"RnaModification",
	"DnaModification",
	"Template",
]