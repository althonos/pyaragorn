import typing
from typing import Literal, Union, Generic
from typing_extensions import Buffer

__version__: str

ARAGORN_VERSION: str
TRANSLATION_TABLES: typing.Set[int]

class Gene:
    @property
    def type(self) -> str: ...
    @property
    def begin(self) -> int: ...
    @property
    def end(self) -> int: ...
    @property
    def length(self) -> int: ...
    @property
    def strand(self) -> Literal[1, -1]: ...
    @property
    def energy(self) -> float: ...
    @property
    def raw_energy(self) -> float: ...
    def sequence(self) -> str: ...
 
class TRNAGene(Gene):
    @property
    def type(self) -> Literal["tRNA"]: ...
    @property
    def amino_acid(self) -> str: ...
    @property
    def amino_acids(self) -> Union[Tuple[str], Tuple[str, str]]: ...
    @property
    def anticodon(self) -> str: ...
    @property
    def anticodon_offset(self) -> int: ...
    @property
    def anticodon_length(self) -> Literal[2, 3, 4]: ...

class TMRNAGene(Gene):
    @property
    def type(self) -> Literal["tmRNA"]: ...
    @property
    def permuted(self) -> bool: ...
    @property
    def orf_offset(self) -> int: ...
    @property
    def orf_length(self) -> int: ...
    def orf(self, include_stop: bool = True) -> str: ...
    def peptide(self, include_stop: bool = True) -> str: ...


class Cursor:
    def __init__(self, obj: Union[str, bytes, bytearray, Buffer]) -> None: ...


G = typing.TypeVar("G", bound=Gene)

class RNAFinder(Generic[G]):
    @typing.overload
    def __init__(
        self: RNAFinder[TRNAGene],
        translation_table: int = 1,
        *,
        tmrna: Literal[False],
        trna: Literal[True] = True,
        linear: bool = False,
        threshold_scale: float = 1.0,
    ) -> None: ...
    @typing.overload
    def __init__(
        self: RNAFinder[TMRNAGene],
        translation_table: int = 1,
        *,
        trna: Literal[False],
        tmrna: Literal[True] = True,
        linear: bool = False,
        threshold_scale: float = 1.0,
    ) -> None: ...
    @typing.overload
    def __init__(
        self: RNAFinder[Union[TRNAGene, TMRNAGene]],
        translation_table: int = 1,
        *,
        trna: bool = True,
        tmrna: bool = True,
        linear: bool = False,
        threshold_scale: float = 1.0,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def translation_table(self) -> int: ...
    @property
    def trna(self) -> bool: ...
    @property
    def tmrna(self) -> bool: ...
    @property
    def linear(self) -> bool: ...
    @property
    def threshold_scale(self) -> float: ...
    @threshold_scale.setter
    def threshold_scale(self, threshold_scale: float) -> None: ...
    def find_rna(
        self,
        sequence: Union[str, bytes, bytearray, Buffer],
    ) -> List[G]: ...
