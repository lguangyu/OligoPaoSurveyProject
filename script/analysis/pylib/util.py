#!/usr/bin/env python3

import dataclasses
import functools
import io


def get_fp(fn_fp, *ka, factory=open, **kw):
	if isinstance(fn_fp, io.IOBase):
		ret = fn_fp
	elif isinstance(fn_fp, str):
		ret = open(fn_fp, *ka, **kw)
	else:
		raise TypeError("fn_fp must be str or io.IOBase, not '%s'"
			% type(fn_fp).__name__)
	return ret


class NonNegFloat(float):
	@functools.wraps(float.__new__)
	def __new__(cls, *ka, **kw):
		new = super().__new__(cls, *ka, **kw)
		if new < 0:
			raise ValueError("%s must be non-negative" % cls.__name__)
		return new


class Fraction(float):
	@functools.wraps(float.__new__)
	def __new__(cls, *ka, **kw):
		new = super().__new__(cls, *ka, **kw)
		if (new < 0) or (new > 1):
			raise ValueError("%s must be between 0 and 1" % cls.__name__)
		return new


class PosInt(int):
	@functools.wraps(int.__new__)
	def __new__(cls, *ka, **kw):
		new = super().__new__(cls, *ka, **kw)
		if new <= 0:
			raise ValueError("%s must be positive" % cls.__name__)
		return new


class NonNegInt(int):
	@functools.wraps(int.__new__)
	def __new__(cls, *ka, **kw):
		new = super().__new__(cls, *ka, **kw)
		if new < 0:
			raise ValueError("%s must be non-negative" % cls.__name__)
		return new


@dataclasses.dataclass
class Taxonomy(object):
	k: str = None  # kingdom
	p: str = None  # phylum
	c: str = None  # class
	o: str = None  # order
	f: str = None  # family
	g: str = None  # genus
	s: str = None  # species

	@classmethod
	def parse(cls, s: str, delimiter=";"):
		levels = ["k", "p", "c", "o", "f", "g", "s"]
		taxa = s.split(delimiter)[:-1]
		new = cls()
		for l, t in zip(levels, taxa):
			setattr(new, l, t)
		return new
