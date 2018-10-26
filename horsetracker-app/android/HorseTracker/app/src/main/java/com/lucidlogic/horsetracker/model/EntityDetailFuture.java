package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;
import android.os.Parcel;
import android.os.Parcelable;

import java.util.Date;
import java.util.List;
import java.util.Map;

/**
 * Created by itowey on 12/08/16.
 */
public class EntityDetailFuture extends BaseObservable implements Parcelable {

    private String id;
    private String sex;
    private String age;
    private String name;
    private String owner;
    private Date collectionDate;
    private Map sire;
    private Map dam;
    private List futEnt;

    public EntityDetailFuture() {
    }

    protected EntityDetailFuture(Parcel in) {
        id = in.readString();
        sex = in.readString();
        age = in.readString();
        name = in.readString();
        owner = in.readString();
        collectionDate = new Date(in.readLong());
        in.readMap(sire, this.getClass().getClassLoader());
        in.readMap(dam, this.getClass().getClassLoader());
        futEnt = in.readArrayList(this.getClass().getClassLoader());
    }

    @Override
    public void writeToParcel(Parcel dest, int flags) {
        dest.writeString(id);
        dest.writeString(sex);
        dest.writeString(age);
        dest.writeString(name);
        dest.writeString(owner);
        dest.writeLong(collectionDate.getTime());
        dest.writeMap(sire);
        dest.writeMap(dam);
        dest.writeList(futEnt);
    }

    @Override
    public int describeContents() {
        return 0;
    }

    public static final Creator<EntityDetailFuture> CREATOR = new Creator<EntityDetailFuture>() {
        @Override
        public EntityDetailFuture createFromParcel(Parcel in) {
            return new EntityDetailFuture(in);
        }

        @Override
        public EntityDetailFuture[] newArray(int size) {
            return new EntityDetailFuture[size];
        }
    };

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getSex() {
        return sex;
    }

    public void setSex(String sex) {
        this.sex = sex;
    }

    public String getAge() {
        return age;
    }

    public void setAge(String age) {
        this.age = age;
    }

    public Map getSire() {
        return sire;
    }

    public void setSire(Map sire) {
        this.sire = sire;
    }

    public Map getDam() {
        return dam;
    }

    public void setDam(Map dam) {
        this.dam = dam;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getOwner() {
        return owner;
    }

    public void setOwner(String owner) {
        this.owner = owner;
    }

    public List getFutEnt() {
        return futEnt;
    }

    public void setFutEnt(List futEnt) {
        this.futEnt = futEnt;
    }

    public Date getCollectionDate() {
        return collectionDate;
    }

    public void setCollectionDate(Date collectionDate) {
        this.collectionDate = collectionDate;
    }
}
