# Reference counting

Reference counting is a memory management technique (pattern), that can help to share objects, using
memory more efficiently and safely. This technique can be applied to objects that are always
transferred via pointers.

## Sharing objects

An object (type) can be shared between many others by having many objects that keep pointers to him.
Sharing an object obviously spares memory and keeps the data always in sync. Often this is an
advantage. This works especially well with immutable objects, but one must be careful with mutable
objects.

The drawback associated to this technique is the deallocation. You need to ensure that you only
deallocate the object if all other pointers to this object are gone to prevent dangling pointers and
you deallocate it before the last pointer to this object is nullified or associated to a different
object to prevent memory leaks.

## Reference counting

Reference counting solves this problem by counting the number of references (pointers) to this
object. At creation, the reference counter is set to 1. Whenever a reference is created to the
object, the counter is increased by 1 and before a reference is removed, the pointer is decreased by
1\. If the reference counter falls to 0, the object is finally released.

This kind of object management creates a contract between callers and callees and is given by the
following three kinds of routines:

- `*_create*`, `*_copy*` routines give back an object with a reference count of one, so that you own
  it and you are responsible of releasing it
- `*_retain` increments the reference count by one, and must be called when you want to add a
  pointer to the object. Now, you are responsible to release the object later.
- `*_release` decrements the reference count by one and deallocates the object if `ref_count` hits
  0\. It must be called when you no longer need your shared copy (i.e. to relinquish ownership).
  Release replaces deallocation.

## Summary

- If you allocate, copy, or retain an object, then you are responsible for releasing the object. If
  you did not allocate, copy, or retain an object, then you should not release it.
- When you receive an object (as the result of a method call), it will normally remain valid until
  the end of your method and the object can be safely returned as a result of your method. If you
  need the object to live longer than this, you must either retain or copy the object.

## Reference counting in CP2k

CP2k employs reference counting for

- objects related to external libraries (BLACS, MPI, ...)
- objects related to meta-data of complex objects (plane wave grids, dense matrices, real-space
  grids, ...)
- memory pools (plane wave grids, dense matrices, ...)
- cell definitions, force environments, wave-function histories, ...

## Sample code

Some sample code to get a feeling of how it works...

```
  ! create the object
  CALL object_create(my_object,...)

  ! Do something with the object
  ...

  ! Reassociate a pointer to an object of the same kind to my_object
  ! Retain my_object first in case of object_pointer has been pointing to my_object
  CALL object_retain(my_object)
  ! Release the object the pointer has been pointing to (release is safe to call with a NULL-pointer)
  CALL object_release(object_pointer)
  object_pointer => my_object

  ! Do something with object and pointer
  ...

  ! If you know the pointer has not been used before, you can just associate it to the object and the order of association and retain does not matter
  NULLIFY(new_object_pointer)
  new_object_pointer => my_object
  CALL object_retain(my_object)

  ! Do something with object and pointers
  ...

  ! Get rid of the object and the pointers in any order you like
  CALL object_release(my_object)
  CALL object_release(object_pointer)
  CALL object_release(new_object_pointer)
```

## Recommendations

- If you implement a new derived type, try to implement the type without reference counting.
  Reference counting should stay the exception and needs to be well justified.
- Apply the correct create (copy)-retain-release-cycle when using a type implementing reference
  counting (indicated by a retain function for the given type).
