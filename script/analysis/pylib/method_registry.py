#!/usr/bin/env python3


class MethodRegistry(dict):
	def __init__(self, *ka, value_type = object):
		if (not isinstance(value_type, object)):
			raise TypeError("value_type must be a class type, not '%s'"\
				% type(value_type).__name__)
		super(MethodRegistry, self).__init__(*ka)
		self.value_type = value_type
		self.default_key = None
		return

	@property
	def default_key(self):
		return self._default_key
	@default_key.setter
	def default_key(self, key: str = None):
		if ((key is not None) and (not isinstance(key, str))):
			raise TypeError("key must be str, not '%s'" % type(key).__name__)
		self._default_key = key
		return

	def register(self, key, *, as_default = False):
		if key in self:
			raise ValueError("key '%s' already exists" % key)
		def decorator(obj):
			if not issubclass(obj, self.value_type):
				raise TypeError("decorated object must be subclass of '%s'"\
					% self.value_type.__name__)
			self[key] = obj
			if as_default:
				self.default_key = key
			return obj
		return decorator

	def get_key_list(self):
		return sorted(self.keys())

	def get(self, key, *ka, **kw):
		cls = self[key]
		return cls(*ka, **kw)
