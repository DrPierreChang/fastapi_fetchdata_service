from abc import ABC, abstractmethod
from typing import Optional, Union

__all__ = (
    "AbstractCache",
    "get_cache",
    "get_access_token_cache",
    "get_refresh_token_cache",
)

from src.core import config


class AbstractCache(ABC):
    def __init__(self, cache_instance):
        self.cache = cache_instance

    @abstractmethod
    def get(self, key: str):
        pass

    @abstractmethod
    def set(
            self,
            key: str,
            value: Union[bytes, str],
            expire: int = config.CACHE_EXPIRE_IN_SECONDS,
    ):
        pass

    @abstractmethod
    def close(self):
        pass

    @abstractmethod
    def delete(self, key: str):
        pass

    @abstractmethod
    def scan_iter(self):
        pass


cache: Optional[AbstractCache] = None
access_token_cache: Optional[AbstractCache] = None
refresh_token_cache: Optional[AbstractCache] = None


# Функция понадобится при внедрении зависимостей
def get_cache() -> AbstractCache:
    return cache


def get_access_token_cache() -> AbstractCache:
    return access_token_cache


def get_refresh_token_cache() -> AbstractCache:
    return refresh_token_cache
