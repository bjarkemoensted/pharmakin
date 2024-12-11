class Registry(dict):
    """Registry class for uniquely registering quantities under a given name.
    Similar to a dict, but throws an error if a name has already been registered."""

    def __setitem__(self, key, value):
        if key in self:
            raise RuntimeError(f"Key {key} is already set (to {value})")
        return super().__setitem__(key, value)
    #


if __name__ == '__main__':
    pass
